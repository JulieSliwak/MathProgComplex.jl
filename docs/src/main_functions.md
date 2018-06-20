# PowSysMod.jl Documentation

## Building one OPF problem for one or several scenarios

Two functions have to be applied.
The first one `load_OPFproblems` converts data contained in `instance_path` in a specific entry `input_type` format into generic data of type OPFProblems.

```@docs
load_OPFproblems(input_type, instance_path::String)
```

OPFProblems is represented thanks to three structures : DataSource, GridStructure, MathematicalProgramming.

The second one `build_globalpb` constructs the OPF problem for all scenarios from generic data `OPFpbs`.

```@docs
build_globalpb!(OPFpbs)
```

A variant: `build_Problem` constructs the OPF problem for a given scenario `scenario` from generic data `OPFpbs`.

```@docs
build_Problem!(OPFpbs, scenario::String)
```

```@repl
include(joinpath(pwd(),"..","..","src_PowSysMod", "PowSysMod_body.jl"))
OPFpbs = load_OPFproblems(MatpowerInput, joinpath(pwd(),"..","..","data_Matpower","matpower","WB2.m"))
OPF = build_Problem!(OPFpbs, "BaseCase")
print(OPF)
OPF = build_globalpb!(OPFpbs)
print(OPF)
```
