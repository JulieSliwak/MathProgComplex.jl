# State of implementation

## The complex hierarchy

The complex hierarchy, working directly on a $POP-\mathbb C$ is not fully implemented and has not been tested. More precisely,

- it is expected that with the current implementation one can get a valid `momentrelaxation::SDPDual{Complex}` and `sosrelaxation::SDPPrimal{Complex}`,
- the structure `SDP_Problem` has to be extended to accept complex values. One way is to make it parametric in `T<:Number`. Care should be taken that the export is also valid for complex coefficients. This should not prove difficult.
- having obtained a `SDP_Problem{Complex}`, it has to be converted to a real SDP instance `SDP_Problem{Float64}`. This conversion, though straightforward in theory, could prove delicate in practice. Hermitian variable matrices $Z_i\in\mathbb H^n$ have to be replaced by real SDP ones $X\in\mathcal S^{2n}, X\succeq0$ with the following additional constraint:

```math
X = \left( \begin{array} \\ A & -B \\ B & A \end{array} \right) \in \mathcal S^{2n}
```

See *Moment/sum-of-squares hierarchy for complex polynomial optimization*, p.3.

## The real hierarchy

### Extensively tested functionalities

Each of these features is individually tested on several instances in the package tests and are therefore fairly reliable:

- Handling of equality constraints,
- Imposing the POP symmetries to set the moments not respecting this symmetry to zero,

### Tested functionalities

Each of these features is individually tested on few instances in the package tests. If behaviour is not inspiring confidence, consider finding more instances and ways to test the feature.

- Clique decomposition: currently tested on the `WB5` $ACOPF$ instance at order 1, `case9` $ACOPF$ instance at orders 1 and 2.
- Handling binary variable: binary variables can be modeled in this continuous framework by a corresponding continuous variable $x^b$ and a constraint $x^b(1-x^b)=0$. This has been added into the module workflow, in such a way that these variables are handled internally and the input $POP$ or clique decomposition is not copied nor modified.

### Expected to be working, but untested functionalities

- Multi-ordered hierarchy: currently working but not tested, by lack of reference data. The only example giving a result on the multi-ordered hierarchy is the `WB5` *complex* SDP relaxation, which should be exact with order 2 relaxation on the power balance and voltage constraints associated with busses $\{4, 5\}$. The example also uses decomposition in two cliques; details are available in *Moment/sum-of-squares hierarchy for complex polynomial optimization*, p.27.

## Working high-level functionalities

- csv export of all information

## Working scripts

- `script_dat_to_hierarchysimple.jl`: For running the hierarchy on OPF problem built from a dat file (see `JulieSliwek/OPFInstances.jl`). Intended to be called as a single process by the batch parallel script `script_dat_to_hierarchyglobal.jl`.
- `script_dat_to_hierarchyglobal.jl`: Runs the hierarchy on several instances, with several settings in parallel using the `qsub` task manager. A pure julia `pmap` implementation could be considered.
- `script_collect_csv.jl`: collect and merge all csv files in the subfolders of the input folder. Serves after running `script_dat_to_hierarchyglobal.jl`.
