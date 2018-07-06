# Benchmark of polynomial module

## `pb_cplx2real` call

Results obtained with `benchamrk/time_cplx2real.jl` code, and the `@btime` macro from `BenchamrkTools` (multiple evaluations).

Improvements consist in removing calls to functions flagged with inefficiency.

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

Use of add! instead of some add in pb_cplx2real (additive algebra).

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

Use of add! in place of all add in pb_cplx2real (additive algebra).

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

### Improvements to import_from_dat

Use of add! in place of all add in import_from_dat (additive algebra).

```bash
Loading module MathProgComplex
elapsed time: 8.838510293 seconds

Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case300.dat")
seek_efficiency() = true
import_from_dat call:
  565.073 ms (2756061 allocations: 120.08 MiB)

pb_cplx2real call:
  2.518 s (20943891 allocations: 930.85 MiB)
```

```bash
Loading module MathProgComplex
elapsed time: 8.5333e-5 seconds

Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case13659pegase.dat")
seek_efficiency() = true
import_from_dat call:
  33.388 s (139179371 allocations: 5.98 GiB)
Base.summarysize(pb_c) = 61674351 #bytes

pb_cplx2real call:
  180.719 s (1032901919 allocations: 44.99 GiB)
Base.summarysize(pb) = 283239022 #bytes
```

### Flag mult. algebra and clean import_from_dat, pb_cplx2real

```bash
Loading module MathProgComplex
elapsed time: 7.0154e-5 seconds

Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case300.dat")
seek_efficiency() = true
import_from_dat call:
  436.200 ms (1841739 allocations: 79.80 MiB)
Base.summarysize(pb_c) = 1382050

pb_cplx2real call:
  2.442 s (20014923 allocations: 889.81 MiB)
Base.summarysize(pb) = 6085556
```
