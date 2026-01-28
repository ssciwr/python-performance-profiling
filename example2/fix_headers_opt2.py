import sys


translation_table = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(translation_table)[::-1]


def reverse_complement_p5(header: str) -> str:
    try:
        h = header.rstrip("\n")
        prefix, last = h.rsplit(":", 1)
        i7, i5 = last.split("+", 1)
        return f"{prefix}:{i7}+{revcomp(i5)}\n"
    except Exception:
        return header  # leave unchanged if malformed


def process_fastq(in_path: str, out_path: str):
    with (
        open(in_path, "r", encoding="ascii") as fin,
        open(out_path, "w", encoding="ascii") as fout,
    ):
        while True:
            header = fin.readline()
            if not header:
                break  # EOF

            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            if header.startswith("@"):
                header = reverse_complement_p5(header)

            fout.write(header)
            fout.write(seq)
            fout.write(plus)
            fout.write(qual)


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <input.fastq> <output.fastq>")

    process_fastq(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
