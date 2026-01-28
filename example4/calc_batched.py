import torch
import torch.nn as nn
from torch.profiler import profile, ProfilerActivity

# Determine device
device = "cuda" if torch.cuda.is_available() else "cpu"
print("device:", device)

# Simple model
model = nn.Sequential(
    nn.Linear(1024, 2048),
    nn.ReLU(),
    nn.Linear(2048, 1024),
).to(device)

# CPU input (kept on CPU so we can force H2D copies each iteration)
x_cpu = torch.randn(512, 1024)

steps = 20

with profile(
    activities=[ProfilerActivity.CPU]
    + ([ProfilerActivity.CUDA] if device == "cuda" else []),
    profile_memory=True,
    record_shapes=True,
    with_stack=True,
) as prof:
    if device == "cuda":
        # CPU -> GPU transfer
        x = x_cpu.to("cuda")
    else:
        x = x_cpu

    for _ in range(steps):
        # GPU compute (or CPU compute if no CUDA)
        y = model(x)

    if device == "cuda":
        # GPU -> CPU transfer
        y_cpu = y.to("cpu")
        _ = float(y_cpu[0, 0])  # touch result so the copy "matters"
        torch.cuda.synchronize()

print(
    prof.key_averages().table(
        sort_by="cuda_time_total" if device == "cuda" else "cpu_time_total"
    )
)
prof.export_chrome_trace("trace.json")
print("Wrote trace.json (open it at https://ui.perfetto.dev)")
