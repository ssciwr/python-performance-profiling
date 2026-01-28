import sys
from Bio.Seq import Seq


def reverse_complement_p5(header):
    try:
        i7, i5 = header.strip().split(":")[-1].split("+")
        i5_rc = str(Seq(i5).reverse_complement())
        return header.replace(i5, i5_rc)
    except Exception:
        return header  # if malformed, skip modification


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
