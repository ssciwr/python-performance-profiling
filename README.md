# [Python Performance Profiling](https://ssciwr.github.io/python-performance-profiling)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Sample Python code to accompany the SSC compact course ["Python Performance Profiling"](https://ssciwr.github.io/python-performance-profiling).

## Getting started

To clone the repo:

```
git clone https://github.com/ssciwr/python-performance-profiling.git
cd python-performance-profiling
```

If using conda, create and activate a new environment:

```
conda create -n python-profiling "python=3.13" -y
conda activate python-profiling
```

Install pytorch for your system following the instructions at https://pytorch.org/get-started/locally/,
for example at time of writing for CUDA 12.8 on linux:

```
pip install torch torchvision
```

Then install the other dependencies:

```
pip install -r requirements.txt
```

## Examples

- [Example 1](example1/README.md): Basic profiling with cProfile, line_profiler.
- [Example 2](example2/README.md): Profiling memory use and performance of a FASTQ file processing script.
- [Example 3](example3/README.md): Profiling a multiprocessing workload.
- [Example 4](example4/README.md): Profiling pytorch GPU code.
