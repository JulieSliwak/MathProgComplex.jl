# Relaxation options

Parameters can be set in the `relaxparams` dictionary of the `RelaxationContext` object. Here is the list of available parameters, along with their default value.

!!! warning

    To be completed, relaxation options will be refactored

## POP parameters

- `:pb_name`: xxxxx
- `:pb_nvar`: xxxxx
- `:pb_nvar_cplx`: xxxxx
- `:pb_nvar_real`: xxxxx
- `:pb_nvar_bin`: xxxxx
- `:pb_nctr`: xxxxx
- `:pb_nctr_eq`: xxxxx
- `:pb_nctr_ineqsimple`: xxxxx
- `:pb_nctr_ineqdouble`: xxxxx
- `:pb_maxpolydeg`: xxxxx
- `:pb_isphaseinv`: xxxxx

## Relaxation build information

- `:slv_mmtrel_t`: xxxxx. Default value is =>-1.0,                    # Moment relaxation time and memory consumption
- `:slv_mmtrel_bytes`:
- `:slv_sosrel_t`: Conversion from moment relaxation to SOS relaxation time (s).
- `:slv_sosrel_bytes`:
- `:slv_fileexport_t`: Export to .sdp files time
- `:slv_fileexport_bytes`:
- `:slv_SDPProblem_t`: Construction of the `SDP_Problem` struct time
- `:slv_SDPProblem_bytes`:

## SDP solver (Mosek) return status and intel

| Option name | Meaning |
| :---------- | :------ |
| `:slv_prosta` | xxxxx. Default value is =>"" |
| `:slv_solsta` | xxxxx. Default value is =>"", |
| `:slv_primobj` | xxxxx. Default value is =>-1.0, |
| `:slv_dualobj` | xxxxx. Default value is =>-1.0, |
| `:slv_primfeas` | xxxxx. Default value is =>-1.0, |
| `:slv_dualfeas` | xxxxx. Default value is =>-1.0, |
| `:slv_solvetime` | xxxxx. Default value is =>-1.0, |

## Options defining the relaxation

- `:opt_hierarchykind`: xxxxx. Default value is =>:Undef,
- `:opt_issparse`: xxxxx. Default value is =>false,                   # Relaxation parameters
- `:opt_multiordered`: xxxxx. Default value is =>false,
- `:opt_globalorder`: xxxxx. Default value is =>-1,
- `:opt_sym_phaseinv`: xxxxx. Default value is =>false,
- `:opt_nb_cliques`: xxxxx. Default value is =>-1,
- `:opt_exportsdp`: xxxxx. Default value is =>0,                      # 0: no, 1: export specified problem to:opt_exportsdppath
- `:opt_exportsdppath`: xxxxx. Default value is =>"SDP_Problem",
- `:opt_msk_maxtime`: xxxxx. Default value is =>-1,                   # Default -1 is no time limit; unit is seconds
- `:opt_outmode`: xxxxx. Default value is =>0,                        # 0: screen, 1: file, 2: both
- `:opt_outlev`: xxxxx. Default value is =>1,                         # 0: none, 1:summary at moment relaxation, sos relaxation, 2: detailled information, 3: full problems
- `:opt_outname`: xxxxx. Default value is =>"momentsos.log",
- `:opt_outcsv`: xxxxx. Default value is =>0,                         # 0: no csv is written, 1: csv is written
- `:opt_outcsvname`: xxxxx. Default value is =>"momentsos_solve.csv")
