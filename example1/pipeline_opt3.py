import random
import time
import re
import unicodedata

import numpy as np

WHITESPACE_RE = re.compile(r"\s+")

# --- High-performance Unicode path (ICU) ---
try:
    import icu  # PyICU

    HAVE_ICU = True
    _NFKD = icu.Normalizer2.getNFKDInstance()

    # Remove combining marks (accents)
    _REMOVE_MARKS = icu.Transliterator.createInstance("[:Nonspacing Mark:] Remove")

    # Remove punctuation (ICU category)
    _REMOVE_PUNCT = icu.Transliterator.createInstance("[:Punctuation:] Remove")

except Exception:
    HAVE_ICU = False
    _NFKD = None
    _REMOVE_MARKS = None
    _REMOVE_PUNCT = None


def generate_records(n):
    records = []
    for i in range(n):
        name = random_name()
        values = generate_values()
        records.append((i, name, values))
    return records


def random_name():
    first = random.choice(["Alice", "Bob", "Carlos", "Diana", "Élodie", "Frédéric"])
    last = random.choice(["Smith", "Müller", "García", "O'Connor", "Núñez"])
    noise = random.choice(["", " ", "  ", "\t"])
    return f"{noise}{first}{noise}{last}{noise}"


def generate_values():
    return [random.random() for _ in range(50)]


def normalize_name(name: str) -> str:
    """
    Normalize a person's name for matching.

    Fast path (if PyICU is available):
    - NFKD normalization (ICU)
    - Case folding (ICU)
    - Remove accents (ICU transliterator)
    - Remove punctuation (ICU transliterator)
    - Normalize whitespace (regex)
    """
    if HAVE_ICU:
        # Unicode normalization (C++ ICU)
        s = _NFKD.normalize(name)

        # Case folding (ICU)
        # foldCase is closer to Python's casefold than lower()
        us = icu.UnicodeString(s)
        us = us.foldCase()
        s = str(us)

        # Remove accents / combining marks (ICU)
        s = _REMOVE_MARKS.transliterate(s)

        # Remove punctuation (ICU)
        s = _REMOVE_PUNCT.transliterate(s)

        # Normalize whitespace
        s = WHITESPACE_RE.sub(" ", s).strip()
        return s

    # --- Fallback: your original unicodedata implementation ---
    name = unicodedata.normalize("NFKD", name)

    # Remove accents
    chars = []
    for c in name:
        if not unicodedata.combining(c):
            chars.append(c)
    name = "".join(chars)

    # Case folding
    name = name.casefold()

    # Remove punctuation
    cleaned = []
    for c in name:
        if c.isalnum() or c.isspace():
            cleaned.append(c)
    name = "".join(cleaned)

    # Normalize whitespace
    name = WHITESPACE_RE.sub(" ", name).strip()

    return name


def mean_and_variance(values):
    n = 0
    mean = 0.0
    M2 = 0.0

    for x in values:
        n += 1
        delta = x - mean
        mean += delta / n
        delta2 = x - mean
        M2 += delta * delta2

    return mean, M2 / n  # population variance


def compute_statistics(records):
    stats = []
    for record_id, name, values in records:
        clean_name = normalize_name(name)
        m, v = mean_and_variance(values)
        stats.append((record_id, clean_name, m, v))
    return stats


def count_similar_record_pairs(stats, mean_threshold=0.01):
    """
    Sort means and count pairs within threshold.
    O(n log n) complexity instead of O(n^2).
    """
    sorted_means = np.sort([s[2] for s in stats])
    count = 0
    for i in range(len(sorted_means)):
        j = np.searchsorted(sorted_means, sorted_means[i] + mean_threshold, side="left")
        count += max(0, j - i - 1)
    return count


def main():
    random.seed(0)

    n = 10000
    records = generate_records(n)

    stats = compute_statistics(records)

    count = count_similar_record_pairs(stats)

    impl = "PyICU" if HAVE_ICU else "unicodedata (fallback)"
    print(f"Normalization impl: {impl}")
    print(f"Found {count} similar record pairs")


if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print(f"Total runtime: {end - start:.2f} seconds")
