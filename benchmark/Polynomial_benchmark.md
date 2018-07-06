# Benchmark of polynomial module

## `pb_cplx2real` call

Results obtained with `benchamrk/time_cplx2real.jl` code, and the `@btime` macro from `BenchamrkTools` (multiple evaluations).

### Start code

```bash
Loading module MathProgComplex
elapsed time: 8.752972617 seconds

Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case300.dat")
import_from_dat call:
  4.381 s (4437915 allocations: 260.50 MiB)

pb_cplx2real call:
  16.043 s (30177358 allocations: 1.52 GiB)
```

### dev_MR version

```bash
Loading module MathProgComplex
elapsed time: 8.822653466 seconds

Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case300.dat")
import_from_dat call:
  4.140 s (4437908 allocations: 261.92 MiB)

pb_cplx2real call:
  15.585 s (29560716 allocations: 1.49 GiB)

pb_cplx2real_add call:
  4.265 s (24470637 allocations: 1.08 GiB)
```

### Further improvements

```bash
Loading module MathProgComplex
elapsed time: 9.677286022 seconds

Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case300.dat")
import_from_dat call:
  4.097 s (4437874 allocations: 261.53 MiB)

pb_cplx2real call:
  16.953 s (29214487 allocations: 1.46 GiB)

pb_cplx2real_add call:
  3.290 s (20943891 allocations: 930.85 MiB)
```