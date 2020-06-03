# Benchmark Report for *PartiallySeparableSolvers*

## Job Properties
* Time of benchmarks:
    - Target: 3 Jun 2020 - 11:08
    - Baseline: 3 Jun 2020 - 11:04
* Package commits:
    - Target: 99afe7
    - Baseline: 99afe7
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

| ID                             | time ratio                   | memory ratio                 |
|--------------------------------|------------------------------|------------------------------|
| `["ros 10 var", "L-BFGS"]`     |                1.07 (5%) :x: |                   1.00 (1%)  |
| `["ros 10 var", "L-SR1"]`      |                1.29 (5%) :x: |                   1.00 (1%)  |
| `["ros 10 var", "P-BFGS"]`     |                1.12 (5%) :x: |                   1.00 (1%)  |
| `["ros 10 var", "P-SR1"]`      | 0.00 (5%) :white_check_mark: | 0.01 (1%) :white_check_mark: |
| `["ros 10 var", "Trunk"]`      |                3.11 (5%) :x: |                   1.00 (1%)  |
| `["ros 10 var", "Trunk_LSR1"]` |                   0.99 (5%)  |                1.36 (1%) :x: |
| `["ros 20 var", "L-BFGS"]`     | 0.69 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 20 var", "L-SR1"]`      | 0.92 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 20 var", "P-BFGS"]`     | 0.82 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 20 var", "P-SR1"]`      | 0.75 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 20 var", "Trunk"]`      | 0.36 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 20 var", "Trunk_LSR1"]` |                   1.00 (5%)  |                1.95 (1%) :x: |
| `["ros 30 var", "L-BFGS"]`     | 0.82 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 30 var", "L-SR1"]`      | 0.92 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 30 var", "P-BFGS"]`     | 0.79 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 30 var", "Trunk"]`      | 0.55 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros 30 var", "Trunk_LSR1"]` |                   1.00 (5%)  |                1.17 (1%) :x: |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["ros 10 var"]`
- `["ros 20 var"]`
- `["ros 30 var"]`

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
       #1  1498 MHz     402453            0       347906      6363328       136656  ticks
       #2  1498 MHz     298640            0       137656      6677171        10015  ticks
       #3  1498 MHz     628546            0       183156      6301765         4750  ticks
       #4  1498 MHz     409296            0       153343      6550828         4031  ticks
       #5  1498 MHz     518250            0       182843      6412375         5578  ticks
       #6  1498 MHz     486218            0       223828      6403421         3468  ticks
       #7  1498 MHz     605187            0       196062      6312218         7265  ticks
       #8  1498 MHz     432093            0       107562      6573796         2328  ticks
       
  Memory: 31.775043487548828 GB (19713.51953125 MB free)
  Uptime: 7113.0 sec
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
       #1  1498 MHz     373703            0       336781      6177515       134078  ticks
       #2  1498 MHz     280406            0       133906      6473468         9937  ticks
       #3  1498 MHz     573468            0       175000      6139312         4531  ticks
       #4  1498 MHz     370359            0       143109      6374296         3562  ticks
       #5  1498 MHz     472312            0       168640      6246828         5187  ticks
       #6  1498 MHz     455890            0       218640      6213234         3343  ticks
       #7  1498 MHz     545984            0       187671      6154109         6968  ticks
       #8  1498 MHz     406875            0       104546      6376328         2312  ticks
       
  Memory: 31.775043487548828 GB (19194.42578125 MB free)
  Uptime: 6887.0 sec
  Load Avg:  0.0  0.0  0.0
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.1 (ORCJIT, cannonlake)
```