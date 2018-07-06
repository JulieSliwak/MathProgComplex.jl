# Benchmark of hierarchy module

Results obtained with `benchamrk/time_hierarchy.jl` code, and the `@btime` macro from `BenchamrkTools` (multiple evaluations).

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
