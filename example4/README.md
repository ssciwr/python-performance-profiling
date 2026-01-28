# Example 4

For this example you need to have PyTorch installed. Follow the instructions at https://pytorch.org/get-started/locally/ to install,
for example at time of writing for CUDA 12.8 on linux:

```
pip install torch torchvision
```

## Running the example

To run the example, execute the following command in your terminal:

```
python calc.py
```

Then view the trace in Chrome by navigating to https://ui.perfetto.dev/ and loading the `result.json` file.

For comparison, the script calc_batch.py does the same calculations but without copying the data to and from the GPU each time.

## Traces

Traces for the two examples are included if you don't have a GPU.
