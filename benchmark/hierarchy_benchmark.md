# Benchmark of hierarchy module

Results obtained with `benchmark/time_hierarchy.jl` code, and the `@btime` macro from `BenchmarkTools` (multiple evaluations).

## Start code

**Warning**: this is a Dict implementation...
This is order 1, no symmetry, on a real problem.

```bash
Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case30pwl.dat")
seek_efficiency() = false
--- time build_sparsity
  3.467 ms (17532 allocations: 889.33 KiB)
--- time build_momentrelaxation
  309.091 ms (2158850 allocations: 100.88 MiB)
--- time build_SOSrelaxation
  16.298 ms (79377 allocations: 2.45 MiB)
--- time read_SDPPrimal
  10.910 ms (156076 allocations: 5.58 MiB)
```

```bash
Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case89pegase.dat")
seek_efficiency() = false
--- time build_sparsity
  22.759 ms (79122 allocations: 3.74 MiB)
--- time build_momentrelaxation
  3.123 s (14651677 allocations: 686.02 MiB)
--- time build_SOSrelaxation
  168.625 ms (632951 allocations: 18.57 MiB)
--- time read_SDPPrimal
  110.278 ms (1212697 allocations: 41.68 MiB)
```

```bash
Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case300.dat")
seek_efficiency() = false
--- time build_sparsity
  99.194 ms (178276 allocations: 8.75 MiB)
--- time build_momentrelaxation
  46.950 s (105821121 allocations: 4.81 GiB)
--- time build_SOSrelaxation
  3.103 s (6725351 allocations: 174.17 MiB)
--- time read_SDPPrimal
  7.919 s (12171381 allocations: 416.78 MiB)
```

## Dict, with unefficient accessor (see `add_index!`)

```bash
Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "WB2.dat")
seek_efficiency() = false
MathProgComplex.DictType = Dict
--- time build_sparsity
  85.743 μs (736 allocations: 40.66 KiB)
--- time build_momentrelaxation
  8.883 ms (54364 allocations: 2.54 MiB)
--- time build_SOSrelaxation
  2.490 ms (2520 allocations: 79.33 KiB)
--- time read_SDPPrimal
  660.100 μs (2794 allocations: 253.88 KiB)


Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case30pwl.dat")
seek_efficiency() = false
MathProgComplex.DictType = Dict
--- time build_sparsity
  3.201 ms (17532 allocations: 889.33 KiB)
--- time build_momentrelaxation
  406.424 ms (2188001 allocations: 101.61 MiB)
--- time build_SOSrelaxation
  98.996 ms (98799 allocations: 2.87 MiB)
--- time read_SDPPrimal
  10.328 ms (156092 allocations: 5.59 MiB)


Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case89pegase.dat")
seek_efficiency() = false
MathProgComplex.DictType = Dict
--- time build_sparsity
  20.176 ms (79122 allocations: 3.74 MiB)
--- time build_momentrelaxation
  3.429 s (14779992 allocations: 689.18 MiB)
--- time build_SOSrelaxation
  889.406 ms (752901 allocations: 21.29 MiB)
--- time read_SDPPrimal
  100.742 ms (1212713 allocations: 41.68 MiB)


Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case300.dat")
seek_efficiency() = false
MathProgComplex.DictType = Dict
--- time build_sparsity
  85.256 ms (178276 allocations: 8.75 MiB)
--- time build_momentrelaxation
  58.382 s (106108603 allocations: 4.81 GiB)
--- time build_SOSrelaxation
  9.706 s (6935010 allocations: 178.74 MiB)
--- time read_SDPPrimal
  7.710 s (12171397 allocations: 416.78 MiB)
```

## SortedDict

```bash
Loading module MathProgComplex
elapsed time: 6.9333e-5 seconds

Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "WB2.dat")
seek_efficiency() = false
MathProgComplex.DictType = DataStructures.SortedDict
--- time build_sparsity
  103.384 μs (1039 allocations: 53.77 KiB)
--- time build_momentrelaxation
  10.272 ms (82467 allocations: 3.81 MiB)
--- time build_SOSrelaxation
  2.368 ms (21277 allocations: 921.58 KiB)
--- time read_SDPPrimal
  696.203 μs (2794 allocations: 253.81 KiB)

Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case30pwl.dat")
seek_efficiency() = false
MathProgComplex.DictType = DataStructures.SortedDict
--- time build_sparsity
  4.355 ms (22949 allocations: 1.10 MiB)
--- time build_momentrelaxation
  1.220 s (8212354 allocations: 376.69 MiB)
--- time build_SOSrelaxation
  476.280 ms (3183845 allocations: 136.97 MiB)
--- time read_SDPPrimal
  11.109 ms (156092 allocations: 5.58 MiB)

Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case89pegase.dat")
seek_efficiency() = false
MathProgComplex.DictType = DataStructures.SortedDict
--- time build_sparsity
  24.062 ms (99058 allocations: 4.58 MiB)
--- time build_momentrelaxation
  18.231 s (80453541 allocations: 3.60 GiB)
--- time build_SOSrelaxation
  6.197 s (29419580 allocations: 1.23 GiB)
--- time read_SDPPrimal
  104.315 ms (1212713 allocations: 41.68 MiB)

Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case300.dat")
seek_efficiency() = false
MathProgComplex.DictType = DataStructures.SortedDict
--- time build_sparsity
  94.083 ms (232356 allocations: 11.06 MiB)
--- time build_momentrelaxation
  249.232 s (988216825 allocations: 44.24 GiB)
--- time build_SOSrelaxation
  66.104 s (287139616 allocations: 12.00 GiB)
--- time read_SDPPrimal
  15.146 s (12171397 allocations: 416.43 MiB)
```

## OrderedDict

```bash
Loading module MathProgComplex
elapsed time: 11.901092386 seconds


Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "WB2.dat")
seek_efficiency() = false
MathProgComplex.DictType = DataStructures.OrderedDict
--- time build_sparsity
  110.359 μs (764 allocations: 40.97 KiB)
--- time build_momentrelaxation
  11.147 ms (56702 allocations: 2.60 MiB)
--- time build_SOSrelaxation
  3.205 ms (3524 allocations: 106.16 KiB)
--- time read_SDPPrimal
  697.434 μs (2794 allocations: 253.88 KiB)


Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case30pwl.dat")
seek_efficiency() = false
MathProgComplex.DictType = DataStructures.OrderedDict
--- time build_sparsity
  4.120 ms (17940 allocations: 898.28 KiB)
--- time build_momentrelaxation
  746.596 ms (2302585 allocations: 104.17 MiB)
--- time build_SOSrelaxation
  135.569 ms (146419 allocations: 4.35 MiB)
--- time read_SDPPrimal
  13.247 ms (156092 allocations: 5.59 MiB)


Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case89pegase.dat")
seek_efficiency() = false
MathProgComplex.DictType = DataStructures.OrderedDict
--- time build_sparsity
  27.488 ms (80329 allocations: 3.75 MiB)
--- time build_momentrelaxation
  6.061 s (15612268 allocations: 706.31 MiB)
--- time build_SOSrelaxation
  1.154 s (1090161 allocations: 31.35 MiB)
--- time read_SDPPrimal
  115.944 ms (1212714 allocations: 41.75 MiB)


Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case300.dat")
seek_efficiency() = false
MathProgComplex.DictType = DataStructures.OrderedDict
--- time build_sparsity
  134.386 ms (182098 allocations: 8.79 MiB)
--- time build_momentrelaxation
  68.519 s (113914430 allocations: 4.96 GiB)
--- time build_SOSrelaxation
  9.550 s (9720947 allocations: 257.92 MiB)
--- time read_SDPPrimal
  6.398 s (12171402 allocations: 417.49 MiB)
```