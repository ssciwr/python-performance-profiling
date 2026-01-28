#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import random
import time
from multiprocessing import Pool
from typing import Iterable, List, Tuple


def work_unit(seed_and_size: Tuple[int, int]) -> Tuple[int, int]:
    """
    One CPU-heavy unit of work.

    Returns:
      (seed, count)
    """
    seed, size = seed_and_size

    # Make each task deterministic but different
    rnd = random.Random(seed)

    # Do some work, where some tasks do more work than others
    start = rnd.randint(10_000, 12_000)
    span = size // 3 + rnd.randint(0, size // 6)
    count = 0
    for n in range(start, start + span):
        r = int(math.isqrt(n))
        for i in range(3, r + 1, 2):
            count += 1

    return seed, count


def run_serial(jobs: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    results = []
    for job in jobs:
        results.append(work_unit(job))
    return results


def run_parallel(
    jobs: Iterable[Tuple[int, int]],
    workers: int,
    chunksize: int = 1,
) -> List[Tuple[int, int]]:
    with Pool(processes=workers) as pool:
        return pool.map(work_unit, jobs, chunksize=chunksize)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Multiprocessing workload for profiling demos."
    )
    parser.add_argument("--size", type=int, default=200_000, help="Work size per task.")
    parser.add_argument("--tasks", type=int, default=50, help="Number of tasks.")
    parser.add_argument(
        "--workers", type=int, default=1, help="Number of worker processes."
    )
    parser.add_argument(
        "--chunksize", type=int, default=1, help="Chunksize for multiprocessing.Pool."
    )
    args = parser.parse_args()

    # Build deterministic jobs list
    jobs = [(i, args.size) for i in range(args.tasks)]

    t0 = time.perf_counter()
    if args.workers > 1:
        run_parallel(jobs, workers=args.workers, chunksize=args.chunksize)
    else:
        run_serial(jobs)
    t1 = time.perf_counter()
    print(f"{args.workers} workers - time: {t1 - t0:.3f}s")


if __name__ == "__main__":
    main()
