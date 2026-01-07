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
conda create -n python-profiling python=3.14 -y
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
