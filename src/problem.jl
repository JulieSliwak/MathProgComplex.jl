export Problem, print

"""
    Problem()

Define an empty mathematical polynomial optimization problem, characterized by
a polynomial `objective`, a dictionnary of `constraints` and a dictionnary of
used variables `variables`.

### Attributes
- `objective::Polynomial` : polynomial criterion.
- `constraints::SortedDict{String, Constraint}` : dictionnary of constraint name to
constraint.
- `variables::SortedDict{String, Type}` : dictionnary of variable name to type.

### Exemple
```julia
julia > x, y, z = ...
julia > pb = Problem()
julia > set_objective(pb, abs2(x+y+z))
julia > add_constraint!(pb, abs2(x) << 1)
julia > add_constraint!(pb, y << 0.5+1im)
julia > add_constraint!(pb, 0-1im << z << 1+0im)
```
"""
mutable struct Problem
  objective::Polynomial
  constraints::SortedDict{String, Constraint}
  variables::SortedDict{String, Type}
end

Problem() = Problem(Polynomial(), SortedDict{String, Constraint}(), SortedDict{String, Type}())


#############################
## Print
#############################
function Base.print(io::IO, pb::Problem)
  print(io, "▶ variables: ")
  for (varName, typ) in pb.variables
    print(io, Variable(varName, typ), " ")
  end
  println(io, "\n▶ objective: ", pb.objective)
  println(io, "▶ constraints: ")
  for (cstrName, cstr) in pb.constraints
    @printf(io, " → %10s: ", cstrName)
    println(io, cstr)
  end
end
