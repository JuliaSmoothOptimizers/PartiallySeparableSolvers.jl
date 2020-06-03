# Benchmark Report for *PartiallySeparableSolvers*

## Job Properties
* Time of benchmarks:
    - Target: 3 Jun 2020 - 10:56
    - Baseline: 3 Jun 2020 - 10:54
* Package commits:
    - Target: 22ef59
    - Baseline: 22ef59
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

| ID                                       | time ratio                   | memory ratio                 |
|------------------------------------------|------------------------------|------------------------------|
| `["ros [10, 20, 30] var", "P-BFGS"]`     | 0.90 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros [10, 20, 30] var", "P-SR1"]`      | 0.85 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros [10, 20, 30] var", "Trunk"]`      | 0.66 (5%) :white_check_mark: |                   1.00 (1%)  |
| `["ros [10, 20, 30] var", "Trunk_LSR1"]` |                   1.00 (5%)  | 0.67 (1%) :white_check_mark: |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["ros [10, 20, 30] var"]`

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
       #1  1498 MHz     323125            0       304234      5775140       121062  ticks
       #2  1498 MHz     243921            0       121109      6037250         9578  ticks
       #3  1498 MHz     485609            0       157265      5759406         4015  ticks
       #4  1498 MHz     309796            0       126609      5965875         3000  ticks
       #5  1498 MHz     400515            0       148843      5852921         4421  ticks
       #6  1498 MHz     404234            0       200812      5797234         2921  ticks
       #7  1498 MHz     456968            0       171953      5773359         6515  ticks
       #8  1498 MHz     354906            0        96765      5950578         2234  ticks
       
  Memory: 31.775043487548828 GB (21379.85546875 MB free)
  Uptime: 6402.0 sec
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
       #1  1498 MHz     309843            0       299453      5690171       120296  ticks
       #2  1498 MHz     238156            0       119265      5941828         9500  ticks
       #3  1498 MHz     457109            0       153734      5688406         3921  ticks
       #4  1498 MHz     294609            0       123546      5881093         2890  ticks
       #5  1498 MHz     379921            0       144593      5774734         4312  ticks
       #6  1498 MHz     385312            0       196515      5717421         2765  ticks
       #7  1498 MHz     431984            0       168015      5699250         6437  ticks
       #8  1498 MHz     339984            0        94781      5864453         2218  ticks
       
  Memory: 31.775043487548828 GB (21968.43359375 MB free)
  Uptime: 6299.0 sec
  Load Avg:  0.0  0.0  0.0
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.1 (ORCJIT, cannonlake)
```