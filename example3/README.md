# Example 3

Make sure you have installed the dependencies by running (in the top level directory of this repository):

```
pip install -r requirements.txt
```

## Running the example

To see the possible arguments for the script, run `python script.py --help`:

```
usage: script.py [-h] [--size SIZE] [--workers WORKERS]

Multiprocessing workload for profiling demos.

options:
  -h, --help            show this help message and exit
  --size SIZE           Work size per task.
  --tasks TASKS         Number of tasks.
  --workers WORKERS     Number of worker processes.
  --chunksize CHUNKSIZE
                        Chunksize for multiprocessing.Pool.
```

To run the script serially (without multiprocessing):

```
python script.py
```

To run the example with multiprocessing enabled (using 4 worker processes):

```
python script.py --workers 4
```

## Profiling

Profile the code using py-spy to make a flamegraph (make sure to include the `-s` flag to capture subprocesses if using multiprocessing):

```
py-spy record -o profile.svg -s -- python script.py --workers=4
```

Use viztracer to generate a trace (you may want to use the `--ignore_c_function` flag to reduce the size of the trace):

```
viztracer --ignore_c_function script.py --workers=4
```

Then view the trace in Chrome by navigating to https://ui.perfetto.dev/ and loading the `result.json` file.

Experiment with changing the size, number of tasks, number of workers, and chunksize to see how it affects the performance.
