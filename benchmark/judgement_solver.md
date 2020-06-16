# Benchmark Report for *PartiallySeparableSolvers*

## Job Properties
* Time of benchmarks:
    - Target: 16 Jun 2020 - 12:11
    - Baseline: 16 Jun 2020 - 11:35
* Package commits:
    - Target: 5f4571
    - Baseline: 5f4571
* Julia commits:
    - Target: 2d5741
    - Baseline: 2d5741
* Julia command flags:
    - Target: None
    - Baseline: None
* Environment variables:
    - Target: None
    - Baseline: None

## Results
A ratio greater than `1.0` denotes a possible regression (marked with :x:), while a ratio less
than `1.0` denotes a possible improvement (marked with :white_check_mark:). Only significant results - results
that indicate possible regressions or improvements - are shown below (thus, an empty table means that all
benchmark results remained invariant between builds).

| ID                               | time ratio                   | memory ratio                 |
|----------------------------------|------------------------------|------------------------------|
| `["ros 100 var", "L-BFGS"]`      | 0.90 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 100 var", "L-SR1"]`       | 0.80 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 100 var", "P-BFGS"]`      | 0.92 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 100 var", "P-BS"]`        | 0.74 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 100 var", "P-SR1"]`       | 0.81 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 100 var", "Trunk"]`       | 0.76 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 1000 var", "L-BFGS"]`     | 0.93 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 1000 var", "L-SR1"]`      | 0.94 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 1000 var", "P-BFGS"]`     | 0.85 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 1000 var", "P-BS"]`       |                1.05 (5%) :x: |                   1.00 (1%)  |
| `["ros 1000 var", "P-SR1"]`      | 0.84 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 1000 var", "Trunk"]`      |                   1.00 (5%)  |                1.08 (1%) :x: |
| `["ros 1000 var", "Trunk_LSR1"]` |                   1.00 (5%)  |                2.53 (1%) :x: |
| `["ros 200 var", "L-BFGS"]`      | 0.88 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 200 var", "L-SR1"]`       | 0.92 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 200 var", "P-BFGS"]`      | 0.92 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 200 var", "P-BS"]`        | 0.95 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 200 var", "P-SR1"]`       | 0.90 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 200 var", "Trunk"]`       | 0.55 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 200 var", "Trunk_LSR1"]`  |                   1.00 (5%)  |                1.56 (1%) :x: |
| `["ros 2000 var", "L-SR1"]`      | 0.71 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 2000 var", "P-BFGS"]`     | 0.93 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 2000 var", "P-BS"]`       | 0.78 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 2000 var", "P-SR1"]`      | 0.72 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 2000 var", "Trunk"]`      |                   1.00 (5%)  |                1.12 (1%) :x: |
| `["ros 2000 var", "Trunk_LSR1"]` |                   1.00 (5%)  | 0.58 (1%) :white_check_mark: |
| `["ros 500 var", "L-BFGS"]`      | 0.37 (5%) :white_check_mark: |                   0.99 (1%)  |
| `["ros 500 var", "L-SR1"]`       | 0.46 (5%) :white_check_mark: | 0.91 (1%) :white_check_mark: |
| `["ros 500 var", "P-BFGS"]`      | 0.80 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 500 var", "P-BS"]`        | 0.88 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 500 var", "P-SR1"]`       | 0.83 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 500 var", "Trunk"]`       |                1.08 (5%) :x: |                   1.00 (1%)  |
| `["ros 500 var", "Trunk_LSR1"]`  |                   1.00 (5%)  | 0.96 (1%) :white_check_mark: |
| `["ros 5000 var", "L-SR1"]`      | 0.94 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 5000 var", "P-BFGS"]`     | 0.93 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 5000 var", "P-BS"]`       |                1.12 (5%) :x: |                   1.00 (1%)  |
| `["ros 5000 var", "Trunk"]`      |                   1.00 (5%)  |                1.10 (1%) :x: |
| `["ros 5000 var", "Trunk_LSR1"]` |                   1.00 (5%)  | 0.99 (1%) :white_check_mark: |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["ros 100 var"]`
- `["ros 1000 var"]`
- `["ros 200 var"]`
- `["ros 2000 var"]`
- `["ros 500 var"]`
- `["ros 5000 var"]`

## Julia versioninfo

### Target
```
Julia Version 1.3.1
Commit 2d5741174c (2019-12-30 21:36 UTC)
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
      Microsoft Windows [version 10.0.18362.476]
  CPU: Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz: 
              speed         user         nice          sys         idle          irq
       #1  1498 MHz    3659609            0      2771703     42726250       996234  ticks
       #2  1498 MHz    4906734            0      1401921     42848937       115687  ticks
       #3  1498 MHz    6629765            0      1617234     40910578        37734  ticks
       #4  1498 MHz    3706250            0      1026343     44424843        31593  ticks
       #5  1498 MHz    4986953            0      1677375     42493109        41046  ticks
       #6  1498 MHz    3321140            0      1153046     44683250        27843  ticks
       #7  1498 MHz    4859812            0      1449562     42848062        35562  ticks
       #8  1498 MHz    4167968            0      1420437     43569015        16906  ticks
       
  Memory: 31.775043487548828 GB (13381.41796875 MB free)
  Uptime: 95268.0 sec
  Load Avg:  0.0  0.0  0.0
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.1 (ORCJIT, cannonlake)
```

### Baseline
```
Julia Version 1.3.1
Commit 2d5741174c (2019-12-30 21:36 UTC)
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
      Microsoft Windows [version 10.0.18362.476]
  CPU: Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz: 
              speed         user         nice          sys         idle          irq
       #1  1498 MHz    3450093            0      2637906     40905062       982578  ticks
       #2  1498 MHz    4754828            0      1344546     40893703       114484  ticks
       #3  1498 MHz    6289078            0      1518171     39185828        36328  ticks
       #4  1498 MHz    3465578            0       954609     42572750        29859  ticks
       #5  1498 MHz    4656625            0      1570609     40765703        39109  ticks
       #6  1498 MHz    2962828            0      1066265     42963828        27375  ticks
       #7  1498 MHz    4373812            0      1333218     41285906        34453  ticks
       #8  1498 MHz    3577156            0      1262515     42153250        16265  ticks
       
  Memory: 31.775043487548828 GB (12301.64453125 MB free)
  Uptime: 93104.0 sec
  Load Avg:  0.0  0.0  0.0
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.1 (ORCJIT, cannonlake)
```