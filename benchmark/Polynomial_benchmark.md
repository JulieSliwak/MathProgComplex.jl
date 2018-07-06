# Benchmark of polynomial module

## `pb_cplx2real` call

Results obtained with `benchamrk/time_cplx2real.jl` code, and the `@btime` macro from `BenchamrkTools` (multiple evaluations).

### Start code

```bash
Working on ("C:\\Users\\gbareilles\\.julia\\v0.6\\OPFInstances\\instances\\data_Matpower\\matpower_QCQP", "case300.dat")
import_from_dat call:
  4.214 s (4437911 allocations: 260.91 MiB)

pb_cplx2real call:
  17.152 s (29560696 allocations: 1.49 GiB)
```
