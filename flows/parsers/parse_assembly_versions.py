"""Daily incremental updates to historical assembly records.

Identifies assembly versions newly superseded since the last run and appends
them to assembly_historical.tsv.  No NCBI fetches are required — data is
copied directly from the previous assembly_current.tsv parse output.

Usage:
    python -m flows.parsers.parse_assembly_versions \\
        --input_path assembly_data_report.jsonl
"""

import csv
import json
import os
from glob import glob
from pathlib import Path
from typing import Optional

from flows.lib.assembly_versions_utils import (
    COL_ASSEMBLY_ID,
    COL_SUPERSEDED_BY,
    COL_SUPERSEDED_BY_VERSION,
    COL_SUPERSEDED_DATE,
    COL_VERSION_STATUS,
    get_accession,
    get_assembly_id,
    open_tsv,
    parse_accession,
    resolve_current_tsv_paths,
    resolve_historical_tsv_path,
)
from flows.lib.conditional_import import flow
from flows.lib.shared_args import INPUT_PATH, YAML_PATH
from flows.lib.shared_args import parse_args as _parse_args
from flows.lib.shared_args import required
from flows.lib.utils import Config, Parser, load_config

DELIMITER = "\t"


def derive_assembly_version_paths(
    input_path: str, config: Optional[Config] = None
) -> tuple[str, str]:
    """Derive previous_tsv and historical_tsv paths from the input JSONL path.

    Both files live alongside the JSONL, following the convention used by
    parse_ncbi_assemblies.  The current-TSV filename — and therefore the name
    of its ".previous" snapshot — comes from config.meta["file_name"] when a
    config is supplied, so nothing here hardcodes it.

    Args:
        input_path (str): Path to the current assembly_data_report.jsonl.
        config (Config, optional): Loaded YAML config naming the current TSV.

    Returns:
        tuple: (previous_tsv, historical_tsv) absolute paths.
    """
    work_dir = os.path.dirname(os.path.abspath(input_path))
    _, previous_tsv, _ = resolve_current_tsv_paths(work_dir, config=config)
    return previous_tsv, resolve_historical_tsv_path(work_dir)


def load_previous_parsed_by_base(previous_tsv: str) -> dict[str, dict[int, dict]]:
    """Load previous parsed results indexed by base accession and version.

    The snapshot is opened through open_tsv, so a gzipped snapshot of a
    gzipped current TSV is read transparently, and the accession is read via
    get_accession rather than assuming a single column name.

    Args:
        previous_tsv (str): Path to assembly_current.tsv from the previous run.

    Returns:
        dict: Nested mapping of base_accession -> version -> row data.
            Returns an empty dict if the file is not found, which is expected
            on the first run after the Phase 0 backfill.
    """
    previous_by_base: dict[str, dict[int, dict]] = {}

    try:
        with open_tsv(previous_tsv) as f:
            for row in csv.DictReader(f, delimiter=DELIMITER):
                accession = get_accession(row)
                if not accession:
                    continue
                base_acc, version = parse_accession(accession)
                if base_acc not in previous_by_base:
                    previous_by_base[base_acc] = {}
                previous_by_base[base_acc][version] = dict(row)
    except FileNotFoundError:
        print(f"Warning: Previous TSV not found: {previous_tsv}")
        print("  This is expected for the first run after the Phase 0 backfill.")
        return {}

    total = sum(len(v) for v in previous_by_base.values())
    print(f"Loaded {total} assemblies from previous parsed results.")
    print(f"  Unique base accessions: {len(previous_by_base)}")

    return previous_by_base


def build_superseded_row(
    previous_row: dict,
    previous_version: int,
    new_accession: str,
    new_version: int,
    release_date: str,
) -> dict:
    """Build a superseded row from a previous row with updated metadata.

    Args:
        previous_row (dict): Row data copied from the previous parsed TSV.
        previous_version (int): Version number of the assembly being superseded.
        new_accession (str): Accession of the assembly that supersedes this one.
        new_version (int): Version number of the superseding assembly.
        release_date (str): Release date of the superseding assembly.

    Emits the canonical camelCase column names declared in
    assembly_historical.types.yaml, so Phase 1 rows are indistinguishable
    from Phase 0 rows in the historical TSV.

    Returns:
        dict: Updated row with versionStatus, assemblyID and the three
            supersession columns set.
    """
    base_acc, _ = parse_accession(get_accession(previous_row))
    row = previous_row.copy()
    row[COL_VERSION_STATUS] = "superseded"
    row[COL_ASSEMBLY_ID] = f"{base_acc}_{previous_version}"
    row[COL_SUPERSEDED_BY] = new_accession
    row[COL_SUPERSEDED_BY_VERSION] = new_version
    row[COL_SUPERSEDED_DATE] = release_date
    return row


def build_missing_version_record(
    base_acc: str,
    missing_version: int,
    new_version: int,
    new_accession: str,
    is_new_series: bool = False,
) -> dict:
    """Build a record describing a version missing from the previous parsed TSV.

    Args:
        base_acc (str): Base accession without version suffix.
        missing_version (int): The version number that could not be found.
        new_version (int): The new version number that triggered this check.
        new_accession (str): Full accession of the new assembly.
        is_new_series (bool): True if this is a new assembly series with no
            prior history at all.

    Returns:
        dict: Record suitable for writing to a missing-versions JSON file.
    """
    record = {
        "base_accession": base_acc,
        "missing_version": missing_version,
        "new_version": new_version,
        "new_accession": new_accession,
    }
    if is_new_series:
        record["note"] = "New assembly series — prior versions may need backfill"
    return record


def identify_newly_superseded(
    new_jsonl: str,
    previous_by_base: dict[str, dict[int, dict]],
) -> tuple[list[dict], list[dict]]:
    """Identify versions that became superseded in the current JSONL update.

    For each assembly with version > 1 in the new JSONL, checks whether the
    immediately preceding version exists in the previous parsed TSV.  Assemblies
    whose predecessor is found are added to the superseded list; those missing a
    predecessor are recorded for optional backfill.

    Args:
        new_jsonl (str): Path to the current assembly_data_report.jsonl.
        previous_by_base (dict): Indexed previous parsed results from
            load_previous_parsed_by_base.

    Returns:
        tuple: (newly_superseded, missing_versions) where each is a list of dicts.
    """
    newly_superseded: list[dict] = []
    missing_versions: list[dict] = []

    with open(new_jsonl) as f:
        for line in f:
            assembly = json.loads(line)
            accession = assembly["accession"]
            base_acc, new_version = parse_accession(accession)

            if new_version <= 1:
                continue

            previous_version = new_version - 1

            if base_acc not in previous_by_base:
                missing_versions.append(build_missing_version_record(
                    base_acc, previous_version, new_version, accession,
                    is_new_series=True,
                ))
                continue

            if previous_version not in previous_by_base[base_acc]:
                missing_versions.append(build_missing_version_record(
                    base_acc, previous_version, new_version, accession,
                ))
                continue

            previous_row = previous_by_base[base_acc][previous_version]
            release_date = assembly.get("releaseDate") or ""
            newly_superseded.append(build_superseded_row(
                previous_row, previous_version, accession, new_version, release_date,
            ))

    return newly_superseded, missing_versions


def merge_fieldnames(
    rows: list[dict], preferred_order: Optional[list[str]] = None
) -> list[str]:
    """Build an output header covering every key present in any row.

    Deriving the header from a single arbitrary row silently drops any column
    that row happens to lack, so the union is taken across all rows.  Columns
    are emitted in preferred_order first (the YAML header order, or the order
    already on disk), then any remaining keys in first-seen order.

    Args:
        rows (list): All row dicts destined for the output file.
        preferred_order (list, optional): Column order to honour where
            possible.

    Returns:
        list: Ordered fieldnames covering the union of all row keys.
    """
    union: dict[str, None] = {}
    for row in rows:
        for key in row:
            union.setdefault(key, None)

    order = list(preferred_order or [])
    ordered = [col for col in order if col in union]
    ordered += [col for col in union if col not in set(order)]
    return ordered


def append_superseded_to_tsv(
    newly_superseded: list[dict],
    historical_tsv: str,
    headers: Optional[list[str]] = None,
) -> None:
    """Merge newly superseded rows into the historical TSV, keyed by assembly ID.

    Reads the existing file if present, merges new rows (new rows take
    precedence on duplicate assembly ID), and writes the combined result back.
    The output header is the union of the columns across all rows, so merging
    Phase 0 and Phase 1 row sets cannot silently drop a column.

    Args:
        newly_superseded (list): Row dicts from identify_newly_superseded.
        historical_tsv (str): Path to assembly_historical.tsv.
        headers (list, optional): Preferred column order, e.g.
            config.headers from assembly_historical.types.yaml.  Defaults to
            the order already present in the file on disk.
    """
    if not newly_superseded:
        print("  No newly superseded versions to add.")
        return

    existing: dict[str, dict] = {}
    file_order: list[str] = []
    historical_path = Path(historical_tsv)

    if historical_path.exists():
        with open_tsv(historical_tsv) as f:
            reader = csv.DictReader(f, delimiter=DELIMITER)
            for row in reader:
                existing[get_assembly_id(row)] = dict(row)
            file_order = list(reader.fieldnames or [])

    for row in newly_superseded:
        existing[get_assembly_id(row)] = row

    fieldnames = merge_fieldnames(
        list(existing.values()), headers or file_order
    )
    with open(historical_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, delimiter=DELIMITER, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(existing.values())

    print(f"  Added {len(newly_superseded)} newly superseded versions.")
    print(f"  Total records in {historical_tsv}: {len(existing)}")


def print_superseded_summary(newly_superseded: list[dict]) -> None:
    """Print a short summary of the newly superseded versions.

    Args:
        newly_superseded (list): Row dicts from identify_newly_superseded.
    """
    if not newly_superseded:
        print("  Found: 0 newly superseded versions.")
        return

    print(f"  Found: {len(newly_superseded)} newly superseded versions.")
    print("  Examples:")
    for row in newly_superseded[:5]:
        print(
            f"    {get_accession(row)} -> "
            f"superseded by v{row[COL_SUPERSEDED_BY_VERSION]}"
        )
    if len(newly_superseded) > 5:
        print(f"    ... and {len(newly_superseded) - 5} more")


def print_missing_versions_warning(missing: list[dict]) -> None:
    """Print a warning listing versions absent from the previous parsed TSV.

    Args:
        missing (list): Missing-version records from identify_newly_superseded.
    """
    if not missing:
        return

    print(f"\n  Warning: {len(missing)} assemblies have missing previous versions.")
    print("  These may need manual backfill:")
    for m in missing[:5]:
        print(
            f"    {m['base_accession']}: "
            f"need v{m['missing_version']}, have v{m['new_version']}"
        )
    if len(missing) > 5:
        print(f"    ... and {len(missing) - 5} more")
    print("\n  To backfill missing versions, run:")
    print("    python -m flows.updaters.update_assembly_versions")


@flow(log_prints=True)
def parse_assembly_versions(
    new_jsonl: str,
    previous_tsv: str,
    historical_tsv: str,
    historical_headers: Optional[list[str]] = None,
) -> dict:
    """Daily incremental update of the historical assembly TSV.

    Called after parse_ncbi_assemblies completes.  Uses the previous parsed
    TSV — no NCBI fetches are made.

    Args:
        new_jsonl (str): Path to the current assembly_data_report.jsonl.
        previous_tsv (str): Path to assembly_current.tsv from the previous run.
        historical_tsv (str): Path to assembly_historical.tsv to update.
        historical_headers (list, optional): Preferred column order for the
            historical TSV.  Defaults to the order already on disk.

    Returns:
        dict: Summary with keys newly_superseded_count, missing_versions_count,
            and missing_versions (list of records).
    """
    separator = "=" * 80
    print(f"\n{separator}")
    print("ASSEMBLY VERSION PARSE")
    print(f"{separator}\n")

    print("[1/3] Loading previous parsed results...")
    previous_by_base = load_previous_parsed_by_base(previous_tsv)

    if not previous_by_base:
        print("  No previous parsed data available — skipping incremental update.")
        print("  This is expected for the first run after the Phase 0 backfill.\n")
        return {
            "newly_superseded_count": 0,
            "missing_versions_count": 0,
            "missing_versions": [],
        }

    print("\n[2/3] Identifying newly superseded versions...")
    newly_superseded, missing = identify_newly_superseded(new_jsonl, previous_by_base)
    print_superseded_summary(newly_superseded)
    print_missing_versions_warning(missing)

    print("\n[3/3] Updating historical TSV...")
    append_superseded_to_tsv(
        newly_superseded, historical_tsv, headers=historical_headers
    )

    print(f"\n{separator}")
    print(
        f"ASSEMBLY VERSION PARSE COMPLETE  "
        f"Superseded: {len(newly_superseded)}  "
        f"Missing: {len(missing)}"
    )
    print(f"{separator}\n")

    return {
        "newly_superseded_count": len(newly_superseded),
        "missing_versions_count": len(missing),
        "missing_versions": missing,
    }


def load_assembly_versions_config(working_yaml: Optional[str]) -> Optional[Config]:
    """Load the parser config, returning None when it is unavailable.

    The config is only used to resolve the current-TSV filename, so a missing
    or unreadable YAML degrades to filename discovery in work_dir rather than
    failing the run.

    Args:
        working_yaml (str, optional): Path to the working YAML file.

    Returns:
        Config or None: Loaded config, or None when it could not be loaded.
    """
    if not working_yaml or not os.path.exists(working_yaml):
        return None
    try:
        return load_config(config_file=working_yaml)
    except Exception as e:
        print(f"  Warning: could not load {working_yaml} ({e});")
        print("  falling back to filename discovery in work_dir.")
        return None


def parse_assembly_versions_wrapper(
    working_yaml: str,
    work_dir: str,
    append: bool,
    data_freeze_path: Optional[str] = None,
    **kwargs,
) -> None:
    """Wrapper matching the fetch_parse_validate parser signature.

    Loads working_yaml so the current-TSV filename — and the ".previous"
    snapshot derived from it — comes from config.meta["file_name"] rather than
    a hardcoded constant, then delegates to parse_assembly_versions.

    Args:
        working_yaml (str): Path to the working YAML file, used to resolve the
            current-TSV filename and its compression.
        work_dir (str): Path to the working directory containing the JSONL,
            the previous TSV, and the historical TSV.
        append (bool): Unused; accepted for pipeline compatibility.
        data_freeze_path (str, optional): Unused; accepted for pipeline
            compatibility.
        **kwargs: Additional keyword arguments passed by the pipeline.
    """
    glob_path = os.path.join(work_dir, "*.jsonl")
    paths = glob(glob_path)
    if not paths:
        raise FileNotFoundError(f"No jsonl file found in {work_dir}")
    if len(paths) > 1:
        raise ValueError(f"More than one jsonl file found in {work_dir}")

    config = load_assembly_versions_config(working_yaml)
    previous_tsv, historical_tsv = derive_assembly_version_paths(
        paths[0], config=config
    )
    results = parse_assembly_versions(
        new_jsonl=paths[0],
        previous_tsv=previous_tsv,
        historical_tsv=historical_tsv,
    )

    if results["missing_versions_count"] > 0:
        missing_json_path = os.path.join(work_dir, "missing_versions.json")
        with open(missing_json_path, "w", encoding="utf-8") as f:
            json.dump(results["missing_versions"], f, indent=2)
        print(f"  Missing versions written to: {missing_json_path}")


def plugin() -> Parser:
    """Register the flow."""
    return Parser(
        name="PARSE_ASSEMBLY_VERSIONS",
        func=parse_assembly_versions_wrapper,
        description="Daily incremental update of historical assembly records.",
    )


if __name__ == "__main__":
    args = _parse_args(
        [required(INPUT_PATH), YAML_PATH],
        description="Daily incremental update of historical assembly records",
    )
    config = load_assembly_versions_config(getattr(args, "yaml_path", None))
    previous_tsv, historical_tsv = derive_assembly_version_paths(
        args.input_path, config=config
    )
    results = parse_assembly_versions(
        new_jsonl=args.input_path,
        previous_tsv=previous_tsv,
        historical_tsv=historical_tsv,
    )
    print(f"Summary: superseded={results['newly_superseded_count']}, "
          f"missing={results['missing_versions_count']}")
    if results["missing_versions_count"] > 0:
        missing_json_path = Path(historical_tsv).parent / "missing_versions.json"
        with open(missing_json_path, "w", encoding="utf-8") as f:
            json.dump(results["missing_versions"], f, indent=2)
        print(
            f"  Action needed: {results['missing_versions_count']} missing versions."
        )
        print(f"  Written to: {missing_json_path}")
        print("  Run: python -m flows.updaters.update_assembly_versions")
