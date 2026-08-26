from types import SimpleNamespace

from flows.parsers.parse_ncbi_assemblies import get_cached_sequence_fields, use_previous_report


def _make_config(previous_row):
    return SimpleNamespace(
        previous_parsed={"GCA_00000001.1": previous_row},
        config={
            "attributes": {
                "total_length": {
                    "header": "totalSequenceLength",
                    "path": "processedAssemblyStats.totalSequenceLength",
                }
            },
            "taxonomy": {
                "genus_taxon_id": {
                    "header": "genusTaxId",
                    "path": "taxonomy.genusTaxId",
                }
            },
        },
    )


def test_release_date_cache_is_normalized_before_reuse():
    processed_report = {
        "processedAssemblyInfo": {"genbankAccession": "GCA_00000001.1"},
        "assemblyInfo": {"releaseDate": "2024-01-01"},
    }
    config = _make_config(
        {
            "genbankAccession": "GCA_00000001.1",
            "releaseDate": " 2024-01-01 ",
            "totalSequenceLength": "12345",
        }
    )

    assert use_previous_report(processed_report, {}, config) is True
    assert get_cached_sequence_fields(processed_report, config) == {"totalSequenceLength": "12345"}


def test_missing_sequence_fields_does_not_count_as_a_valid_cache_hit():
    processed_report = {
        "processedAssemblyInfo": {"genbankAccession": "GCA_00000001.1"},
        "assemblyInfo": {"releaseDate": "2024-01-01"},
    }
    config = _make_config(
        {
            "genbankAccession": "GCA_00000001.1",
            "releaseDate": "2024-01-01",
            "genusTaxId": "9606",
        }
    )

    assert use_previous_report(processed_report, {}, config) is False
    assert get_cached_sequence_fields(processed_report, config) is None
