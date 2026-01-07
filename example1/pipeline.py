import random
import time
import unicodedata
import re


WHITESPACE_RE = re.compile(r"\s+")


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


def normalize_name(name):
    """
    Normalize a person's name for matching.

    Steps:
    - Unicode normalization
    - Case folding
    - Remove accents
    - Remove punctuation
    - Normalize whitespace
    """
    # Unicode normalization
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


def mean(values):
    total = 0.0
    for v in values:
        total += v
    return total / len(values)


def variance(values):
    m = mean(values)
    total = 0.0
    for v in values:
        diff = v - m
        total += diff * diff
    return total / len(values)


def compute_statistics(records):
    stats = []
    for record_id, name, values in records:
        clean_name = normalize_name(name)
        m = mean(values)
        v = variance(values)
        stats.append((record_id, clean_name, m, v))
    return stats


def count_similar_record_pairs(stats, mean_threshold=0.01):
    """
    Naïve similarity search: count record pairs whose means are close.
    """
    count = 0
    for i in range(len(stats)):
        id1, name1, mean1, var1 = stats[i]
        for j in range(i + 1, len(stats)):
            id2, name2, mean2, var2 = stats[j]
            if abs(mean1 - mean2) < mean_threshold:
                count += 1
    return count


def main():
    random.seed(0)

    n = 10000
    records = generate_records(n)

    stats = compute_statistics(records)

    count = count_similar_record_pairs(stats)

    print(f"Found {count} similar record pairs")


if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print(f"Total runtime: {end - start:.2f} seconds")
