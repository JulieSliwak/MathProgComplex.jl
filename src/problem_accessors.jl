###############################################################################
#### Variables
###############################################################################
get_variables(pb::Problem) = pb.variables

function get_variabletype(pb::Problem, varName::String)
  if !haskey(pb.variables, varName)
    error("get_variable(): Current Problem has no variable named ", varName)
  end
  pb.variables[varName]
end

has_variable(pb::Problem, var::Variable) = haskey(pb.variables, var.name) && pb.variables[var.name] == var.kind


function add_variable!(pb::Problem, x::Pair{String, T}) where T
  if haskey(pb.variables, x[1]) && !(x[2] <: pb.variables[x[1]])
    error("add_variable!(): Attempting to add variable ", x.name, " (", x.kind, ") when ", x, " (", get_variabletype(pb, x.name), ") already exists.")
  end
  pb.variables[x[1]] = x[2]
end

"""
  add_variable!(pb, x::Variable)

  Add the `x` variable to `pb` set of variables.
"""
function add_variable!(pb::Problem, x::Variable)
  add_variable!(pb, x.name => x.kind)
end

"""
  add_variables!(pb, expo::Exponent)

  Add variables from `expo` exponent to `pb` set of variables.
"""
function add_variables!(pb::Problem, expo::Exponent)
  for var in keys(expo)
    add_variable!(pb, var)
  end
end

"""
  add_variables!(pb, p::Polynomial)

  Add variables from `p` polynomial to `pb` set of variables.
"""
function add_variables!(pb::Problem, p::Polynomial)
  vars = Set{Variable}()
  for expo in keys(p)
    for var in keys(expo)
      push!(vars, var)
    end
  end
  for var in vars
    add_variable!(pb, var)
  end
  return
end

"""
  add_variables!(pb, ctr::Constraint)

  Add variables from `ctr` constraint to `pb` set of variables.
"""
function add_variables!(pb::Problem, ctr::Constraint)
  add_variables!(pb, ctr.p)
  return
end

###############################################################################
#### Objective
###############################################################################
get_objective(pb::Problem) = pb.objective

"""
get_objective(pb::Problem, pt::Point)

Return the polynomial objective evaluated at `pt`.
"""
get_objective(pb::Problem, pt::Point) = evaluate(pb.objective, pt)

function set_objective!(pb::Problem, p::Polynomial)
  add_variables!(pb, p)
  pb.objective = p
  return
end


###############################################################################
#### Constraints
###############################################################################
has_constraint(pb::Problem, cstrName::String) = haskey(pb.constraints, cstrName)

get_constraints(pb::Problem) = pb.constraints

"""
  cstr = get_constraint(pb::Problem, cstrname::String)

Return the `cstrname` constraint from `pb`.
"""
get_constraint(pb::Problem, cstrname::String) = pb.constraints[cstrname]

"""
  cstr = get_constraint(pb::Problem, pt::Point)

Return the point (dict) of constraints names to their body's value when
evaluated at `pt`.
"""
function get_constraints(pb::Problem, pt::Point)
  img = Point()
  for (cstrName, cstr) in pb.constraints
    img[Variable(cstrName, Complex)] = evaluate(cstr.p, pt)
  end
  return img
end

"""
  add_constraint!(pb::Problem, cstrName::String, cstr::Constraint)

Add the constraint `cstr` under the `cstrname` name in the `pb` problem.
"""
function add_constraint!(pb::Problem, cstrName::String, cstr::Constraint)
  if haskey(pb.constraints, cstrName)
    warn("add_constraint!(): A constraint with that name already exists ($cstrName)")
  end
  if real(cstr.ub) < real(cstr.lb) || imag(cstr.ub) < imag(cstr.lb)
    warn("add_constraint!(): ", cstrName, " Lower bound is higher than upper bound ($(cstr.lb) - $(cstr.ub))")
  end
  add_variables!(pb, cstr)
  pb.constraints[cstrName] = cstr
end

rm_constraint!(pb::Problem, cstrName::String) = pop!(pb.constraints, cstrName)

"""
  pt = get_slacks(pb::Problem, pt::Point)

Return a point associating a constraint name to its slack, i.e. the minimum
algebraic distance between the real and imaginary parts of the body's value at
`pt` and its bounds.
"""
function get_slacks(pb::Problem, pt::Point)
  var_arr = Variable[]
  val_arr = Complex[]
  for (cstrName, cstr) in pb.constraints
    val = evaluate(cstr.p, pt)
    isa(val, Number) || error("get_slacks(): constraint $cstrName not fully evaluated at provided point.\nEvaluated value is $val.")
    push!(var_arr, Variable(cstrName, Complex))
    push!(val_arr, min(real(val-cstr.lb), real(cstr.ub-val)) + min(imag(val-cstr.lb), imag(cstr.ub-val))*im)
  end
  return Point(var_arr, val_arr)
end

function get_relative_slacks(pb::Problem, pt::Point)
  var_arr = Variable[]
  val_arr = Complex[]
  for (cstrName, cstr) in pb.constraints
    val = evaluate(cstr.p, pt)
    lb = cstr.lb
    ub = cstr.ub
    isa(val, Number) || error("get_slacks(): constraint $cstrName not fully evaluated at provided point.\nEvaluated value is $val.")
    push!(var_arr, Variable(cstrName, Complex))
    if lb==ub
      infeas = abs(real(val-ub))/max(abs(real(ub)),1) + im *abs(imag(val-ub))/max(abs(imag(ub)),1)
    else
      infeas_ub = max(0,real(val-ub))/max(abs(real(ub)),1) + im *max(0,imag(val-ub))/max(abs(imag(ub)),1)
      infeas_lb = max(0,real(lb-val))/max(abs(real(ub)),1) + im *max(0,imag(lb-val))/max(abs(imag(ub)),1)
      infeas = max(real(infeas_ub), real(infeas_lb)) + im * max(imag(infeas_ub), imag(infeas_lb))
    end
    push!(val_arr, infeas)
  end
  return Point(var_arr, val_arr)
end


function get_minslack(pb::Problem, pt::Point)
  minSlack = +Inf
  minCstrName = ""
  slacks = get_slacks(pb, pt)
  for (cstrName, slack) in slacks
    if minSlack > min(real(slack), imag(slack))
      minSlack = min(real(slack), imag(slack))
      minCstrName = cstrName
    end
  end
  return minSlack, minCstrName
end

function get_relativemaxslack(pb::Problem, pt::Point)
  maxRelSlack = - Inf
  maxCstrName = ""
  slacks = get_relative_slacks(pb, pt)
  for (cstrName, slack) in slacks
    max_slack = max(real(slack), imag(slack))
    if maxRelSlack < max_slack
      maxRelSlack = max_slack
      maxCstrName = cstrName
    end
  end
  return maxRelSlack, maxCstrName
end


function get_nb_infeasible_ctr_by_ctrtype(pb::Problem, pt::Point, epsilon::Float64)
  slacks = get_relative_slacks(pb, pt)
  nb_infeas_per_ctrtype = Dict{Tuple{String,String},Int64}()
  for (cstrName, slack) in slacks
    if real(slack) > epsilon || imag(slack) > epsilon
      ctr_elemid = split(cstrName.name, '_')
      scenario = ctr_elemid[1]
      bus_or_link = ctr_elemid[2]
      im_or_re = String(ctr_elemid[end])
      ctrtype = String(ctr_elemid[end-1])
      if ismatch(r"BALANCE", ctrtype)
        ctrtype = "BALANCE"
      end
      if !haskey(nb_infeas_per_ctrtype, (ctrtype,im_or_re))
        nb_infeas_per_ctrtype[(ctrtype,im_or_re)] = 0
      end
      nb_infeas_per_ctrtype[(ctrtype, im_or_re)] +=1
    end
  end
  return nb_infeas_per_ctrtype
end
