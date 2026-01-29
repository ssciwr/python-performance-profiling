# Example 2

Make sure you have installed the dependencies by running (in the top level directory of this repository):

```
pip install -r requirements.txt
```

## Running the example

To run the example, first generate a test.fastq example data file by executing:

```
bash generate_data.sh
```

Then you can run the fix_headers script with this file as input:

```
python fix_headers.py test.fastq out.fastq
```

It will read the input FASTQ file, fix the headers, and write the output to `out.fastq`.

## Profiling

Profile the code in the same way as for [example 1](../example1/README.md).

Also try using pyinstrument and py-spy, and compare them with cProfile.

## memray

To profile memory usage with memray, run:

```
memray run -o memray.bin pipeline.py
```
You can then visualize the memory profiling results using memray's built-in viewer:

```
memray flamegraph memray.bin
```


## Benchmarks

Using [hyperfine](https://github.com/sharkdp/hyperfine) to benchmark the scripts,
the first version takes about half a second to process the test data:

```
Benchmark 1: python fix_headers.py test.fastq out.fastq
  Time (mean ± σ):     521.6 ms ±  18.3 ms    [User: 437.0 ms, System: 80.7 ms]
  Range (min … max):   508.2 ms … 566.2 ms    10 runs
```

and the optimised version is about twice as fast:

```
Benchmark 1: python fix_headers_opt2.py test.fastq out.fastq
  Time (mean ± σ):     252.1 ms ±   5.9 ms    [User: 203.7 ms, System: 46.9 ms]
  Range (min … max):   247.7 ms … 269.5 ms    12 runs
```

It would probably be difficult to optimise this much more in pure Python - the next step would
be to implement this in a compiled language like C or Rust, as is done in the tool [fastq-fix-i5](https://crates.io/crates/fastq-fix-i5),
which is about 5 times faster than our optimised Python script (and most of that time is spent in file I/O):

```
Benchmark 1: fastq-fix-i5 < test.fastq > out.fastq
  Time (mean ± σ):      53.0 ms ±   6.5 ms    [User: 15.3 ms, System: 35.0 ms]
  Range (min … max):    45.7 ms …  86.8 ms    49 runs
```
