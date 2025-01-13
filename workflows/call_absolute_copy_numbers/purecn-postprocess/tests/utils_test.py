import pandas as pd
from pytest import approx

from purecn_postprocess.types import Loh
from purecn_postprocess.utils import calculate_cin, call_wgd, type_data_frame


class TestCalculateCin:
    def test_allele_nonspecific_normal(self):
        loh = type_data_frame(
            pd.DataFrame(
                [
                    # ref
                    {"size": 300, "c": 2, "m": 1, "type": None},
                    {"size": 500, "c": 2, "m": 1, "type": None},
                    # aberrant
                    {"size": 200, "c": 3, "m": 1, "type": None},
                ]
            ),
            Loh,
        )

        observed = calculate_cin(loh, allele_specific=False, reference_state="normal")
        expected = 0.2
        assert observed == approx(expected)

    def test_allele_specific_normal(self):
        loh = type_data_frame(
            pd.DataFrame(
                [
                    # ref
                    {"size": 300, "c": 2, "m": 1, "type": None},
                    {"size": 500, "c": 2, "m": 1, "type": None},
                    # aberrant
                    {"size": 100, "c": 3, "m": 1, "type": None},
                    {"size": 100, "c": 2, "m": 0, "type": None},
                    # dropped
                    {"size": 100, "c": 4, "m": None, "type": None},
                ]
            ),
            Loh,
        )

        observed = calculate_cin(loh, allele_specific=True, reference_state="normal")
        expected = 0.2
        assert observed == approx(expected)

    def test_allele_nonspecific_dominant(self):
        loh = type_data_frame(
            pd.DataFrame(
                [
                    # ref
                    {"size": 300, "c": 5, "m": 0, "type": None},
                    {"size": 500, "c": 5, "m": 1, "type": None},
                    # aberrant
                    {"size": 200, "c": 3, "m": 3, "type": None},
                ]
            ),
            Loh,
        )

        observed = calculate_cin(loh, allele_specific=False, reference_state="dominant")
        expected = 0.2
        assert observed == approx(expected)

    def test_allele_specific_dominant(self):
        loh = type_data_frame(
            pd.DataFrame(
                [
                    # ref
                    {"size": 300, "c": 5, "m": 1, "type": None},
                    {"size": 500, "c": 5, "m": 1, "type": None},
                    # aberrant
                    {"size": 100, "c": 3, "m": 1, "type": None},
                    {"size": 100, "c": 2, "m": 0, "type": None},
                    # dropped
                    {"size": 100, "c": 5, "m": None, "type": None},
                ]
            ),
            Loh,
        )

        observed = calculate_cin(loh, allele_specific=True, reference_state="dominant")
        expected = 0.2
        assert observed == approx(expected)


class TestCallWgd:
    def test_not_wgd(self):
        loh = type_data_frame(
            pd.DataFrame(
                [
                    # support
                    {"size": 100, "c": 2, "m": 1, "type": "LOH"},
                    {"size": 200, "c": 2, "m": 1, "type": "COPY-NEUTRAL LOH"},
                    {"size": 200, "c": 2, "m": 1, "type": "WHOLE ARM COPY-NEUTRAL LOH"},
                    # reject
                    {"size": 300, "c": 2, "m": 1, "type": "other"},
                    {"size": 200, "c": 2, "m": 1, "type": None},
                ]
            ),
            Loh,
        )

        assert not call_wgd(loh, ploidy=1.9)

    def test_wgd(self):
        loh = type_data_frame(
            pd.DataFrame(
                [
                    # support
                    {"size": 100, "c": 2, "m": 1, "type": "LOH"},
                    {"size": 200, "c": 2, "m": 1, "type": "COPY-NEUTRAL LOH"},
                    {"size": 200, "c": 2, "m": 1, "type": "WHOLE ARM COPY-NEUTRAL LOH"},
                    # reject
                    {"size": 1000, "c": 2, "m": 1, "type": "other"},
                    {"size": 1000, "c": 2, "m": 1, "type": None},
                ]
            ),
            Loh,
        )

        assert not call_wgd(loh, ploidy=2.0)
