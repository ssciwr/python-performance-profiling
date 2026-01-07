import pytest

from pipeline import (
    normalize_name,
    mean,
    variance,
    compute_statistics,
    count_similar_record_pairs,
)


def test_normalize_name_basic():
    assert normalize_name(" Alice   Smith ") == "alice smith"


def test_normalize_name_accents_and_case():
    assert normalize_name("Élodie Müller") == "elodie muller"


def test_normalize_name_punctuation():
    assert normalize_name("O'Connor!!!") == "oconnor"


def test_mean_simple():
    values = [1.0, 2.0, 3.0]
    assert mean(values) == pytest.approx(2.0)


def test_variance_simple():
    values = [1.0, 2.0, 3.0]
    # Population variance
    expected = 2.0 / 3.0
    assert variance(values) == pytest.approx(expected)


def test_compute_statistics_structure():
    records = [
        (1, " Alice Smith ", [1.0, 2.0, 3.0]),
        (2, "Bob Smith", [2.0, 3.0, 4.0]),
    ]

    stats = compute_statistics(records)

    assert len(stats) == 2

    record_id, name, mean_value, var_value = stats[0]

    assert record_id == 1
    assert name == "alice smith"
    assert isinstance(mean_value, float)
    assert isinstance(var_value, float)


def test_count_similar_record_pairs():
    stats = [
        (1, "alice smith", 1.0, 0.1),
        (2, "bob smith", 1.005, 0.2),
        (3, "carol smith", 2.0, 0.3),
        (4, "jim jones", 1.19, 0.6),
    ]

    assert count_similar_record_pairs(stats, mean_threshold=0.00) == 0  # no pairs
    assert count_similar_record_pairs(stats, mean_threshold=0.01) == 1  # 1-2
    assert count_similar_record_pairs(stats, mean_threshold=0.189) == 2  # 1-2, 2-4
    assert count_similar_record_pairs(stats, mean_threshold=0.191) == 3  # 1-2, 2-4, 1-4
    assert count_similar_record_pairs(stats, mean_threshold=10) == 6  # all pairs


def test_test_count_similar_record_pairs_empty():
    assert count_similar_record_pairs([]) == 0


def test_test_count_similar_record_pairs_single():
    stats = [(1, "alice", 1.0, 0.1)]
    assert count_similar_record_pairs(stats) == 0
