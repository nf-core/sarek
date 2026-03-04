# Benchmark Results

## Test Conditions

| Property | Value |
|----------|-------|
| CPU | {model}, {cores}C/{threads}T |
| RAM | {size} |
| OS | {distro} {version} |
| Storage | {type} |
| Rewritten tool | {name} v{version} (commit {hash}) |
| Upstream tool(s) | {name} v{version} |
| Date | {YYYY-MM-DD} |

## Input Data

| Dataset | Description | Size | Records |
|---------|-------------|------|---------|
| Small | {description} | {file size} | {read/record count} |
| Large | {description} | {file size} | {read/record count} |

## Results

### Wall Clock Time

| Input | Tool | Threads | Time (mean +/- std) | Speedup |
|-------|------|---------|---------------------|---------|
| Small | Upstream | 1 | {time} | 1x |
| Small | Rewritten | 1 | {time} | {N}x |
| Large | Upstream | 1 | {time} | 1x |
| Large | Rewritten | 1 | {time} | {N}x |
| Large | Rewritten | 4 | {time} | {N}x |
| Large | Rewritten | 8 | {time} | {N}x |

### Peak Memory (RSS)

| Input | Tool | Peak RSS |
|-------|------|----------|
| Small | Upstream | {size} |
| Small | Rewritten | {size} |
| Large | Upstream | {size} |
| Large | Rewritten | {size} |

### Thread Scaling (Large Input)

| Threads | Time | Speedup vs 1-thread | Efficiency |
|---------|------|----------------------|------------|
| 1 | {time} | 1.0x | 100% |
| 2 | {time} | {N}x | {N}% |
| 4 | {time} | {N}x | {N}% |
| 8 | {time} | {N}x | {N}% |

## Methodology

- Wall clock times measured with `hyperfine` (3-5 runs, 1 warmup)
- Peak RSS from `/usr/bin/time -v`
- System was idle during benchmarking
- Speedup = upstream_time / rewritten_time
- Thread efficiency = (speedup / num_threads) * 100%

## Notes

{Any relevant observations about the benchmark results}
