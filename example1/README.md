# Example 1

Make sure you have installed the dependencies by running (in the top level directory of this repository):

```
pip install -r requirements.txt
```

## Running the example

To run the example, execute the following command in your terminal:

```
python pipeline.py
```

It does some (made up) data processing, and should take around 10 seconds to run.
You can change the value of `n` in the code to make it run for longer or shorter.
You can also run the test suite if you like:

```
pytest
```

## cProfile

To profile the example using cProfile, run:

```
python -m cProfile -o profile.prof pipeline.py
```

You can then visualize the profiling results using SnakeViz:

```
snakeviz profile.prof
```
This will open a web browser with an interactive visualization of the profiling data.

Or you can generate a call graph using gprof2dot and Graphviz:

```
gprof2dot -f pstats profile.prof | dot -Tsvg -o callgraph.svg
```
This will create a file `callgraph.svg` containing the call graph visualization.

## line_profiler

To use line_profiler, add `from line_profiler import profile` to the imports in the script, then
add the `@profile` decorator to the functions you want to profile in `pipeline.py`.

Then run the script with the `LINE_PROFILE=1` env var set:

```
LINE_PROFILE=1 python pipeline.py
```
