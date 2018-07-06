export cplx2real, pb_cplx2real, real2cplx
export pb_cplx2real_add


varname_cplx2real(varname::String) = (varname*"_Re",varname*"_Im")


"""
	pb = pb_cplx2real(pb_C::Problem)

Convert a complex Polynomial Optimization Problem into a real polynomial
problem in real variables.
"""
function pb_cplx2real(pb_C::Problem)
	pb = Problem()
	for (varName, varType) in get_variables(pb_C)
		if varType <: Complex
			varName_real, varName_imag = varname_cplx2real(varName)
			add_variable!(pb, Variable(varName_real, Real))
			add_variable!(pb, Variable(varName_imag, Real))
		else
			add_variable!(pb, Variable(varName, varType))
		end
	end

	realPart, imagPart = cplx2real(pb_C.objective)
	set_objective!(pb, realPart)

	for (cstrName, cstr) in get_constraints(pb_C)
		realPart, imagPart = cplx2real(cstr.p)
		cstrName_real, cstrName_imag = varname_cplx2real(cstrName)


		if length(realPart) != 0
			# Set bounds to proper infty, easier to detect... TODO : Missing attribute ctr_kind ?
			lb = real(cstr.lb)==-Inf ? -Inf-im*Inf : real(cstr.lb)
			ub = real(cstr.ub)== Inf ? +Inf+im*Inf : real(cstr.ub)
			cstrreal = lb << realPart << ub
			cstr.precond != :none && (cstrreal.precond = cstr.precond)
			add_constraint!(pb, cstrName_real, cstrreal)
		end
		if length(imagPart) != 0
			lb = imag(cstr.lb)==-Inf ? -Inf-im*Inf : imag(cstr.lb)
			ub = imag(cstr.ub)== Inf ? +Inf+im*Inf : imag(cstr.ub)
			cstrimag = lb << imagPart << ub
			cstr.precond != :none && (cstrimag.precond = cstr.precond)
			add_constraint!(pb, cstrName_imag, cstrimag)
		end
	end
	return pb
end


function pb_cplx2real_add(pb_C::Problem)
  # println("---- pb_cplx2real_add ----")
  pb = Problem()
  # println("variables")
  # tic()
  for (varName, varType) in get_variables(pb_C)
    if varType <: Complex
      varName_real, varName_imag = varname_cplx2real(varName)
      add_variable!(pb, Variable(varName_real, Real))
      add_variable!(pb, Variable(varName_imag, Real))
    else
      add_variable!(pb, Variable(varName, varType))
    end
  end
  # toc()
  # println("objective")
  # tic()
  (realPart, imagPart) = cplx2real_add(pb_C.objective)
  set_objective!(pb, realPart)
  # toc()
  # println("constraints")
  # println("nb = ", length(get_constraints(pb_C)))
  # tic()
  constraints = get_constraints(pb_C)
  for (cstrName, cstr) in constraints
    realPart, imagPart = cplx2real_add(cstr.p)
    cstrName_real, cstrName_imag = varname_cplx2real(cstrName)


    if length(realPart) != 0
      # Set bounds to proper infty, easier to detect... TODO : Missing attribute ctr_kind ?
      lb = real(cstr.lb)==-Inf ? -Inf-im*Inf : real(cstr.lb)
      ub = real(cstr.ub)== Inf ? +Inf+im*Inf : real(cstr.ub)
      cstrreal = lb << realPart << ub
      cstr.precond != :none && (cstrreal.precond = cstr.precond)
      add_constraint!(pb, cstrName_real, cstrreal)
    end
    if length(imagPart) != 0
      lb = imag(cstr.lb)==-Inf ? -Inf-im*Inf : imag(cstr.lb)
      ub = imag(cstr.ub)== Inf ? +Inf+im*Inf : imag(cstr.ub)
      cstrimag = lb << imagPart << ub
      cstr.precond != :none && (cstrimag.precond = cstr.precond)
      add_constraint!(pb, cstrName_imag, cstrimag)
    end
  end
  # toc()
  # println("--------------------")
  return pb
end

# Conversion of all complex variables to real ones in the Poly structures

"""
	realPart, imagPart = cplx2real(expo::Exponent)

	Convert a complex Exponent in complex variables into `realPart` and
	`imagPart` polynomials of twice as many variables, real and imag parts of
	`expo` variables. Done recursively with the `cplx2real_rec` function.
"""
function cplx2real(expo::Exponent)
	vars, inds = collect(keys(expo.expo)), collect(values(expo.expo))
	return cplx2real_rec(vars, inds, Polynomial()+1, Polynomial()+0, length(expo)+1, Degree(0,0))
end

"""
realPart, imagPart = cplx2real_rec(vars::Array{Variable}, degs::Array{Degree}, realPart::Polynomial, imagPart::Polynomial, cur_ind::Int, cur_deg::Degree)

Transform recursively the complex exponent represented by the `vars` and `degs`
arrays into its real and imag parts, functions of its imag and real part
variables.
`cur_ind` decreases to 0, `cur_deg` decreases to Degree(0,0) for each step of
cur_ind. Terminaison case is reached at 0, Degree(0,0).
Initial arrays `vars` and `degs` are read only.

### Arguments
- vars::Array{Variable}
- degs::Array{Degree}
- realPart::Polynomial
- imagPart::Polynomial
- cur_ind::Int
- cur_deg::Degree
"""
function cplx2real_rec(vars::Array{Variable}, degs::Array{Degree}, realPart::Polynomial, imagPart::Polynomial, cur_ind::Int, cur_deg::Degree)
	## Final case:
	if cur_ind == 1 && cur_deg == Degree(0,0)
		return (realPart, imagPart)
	## One less variable to deal with:
	elseif cur_deg == Degree(0,0)
		return cplx2real_rec(vars, degs, realPart, imagPart, cur_ind-1, degs[cur_ind-1])
	## Recursion rule, decrease the current variable exponent until it reaches Degree(0,0):
	else
		var = vars[cur_ind]
		if iscomplex(var)
			var_real, var_imag = varname_cplx2real(var.name)
			var_R, var_I = Variable(var_real, Real), Variable(var_imag, Real)
			if cur_deg.explvar > 0
				return cplx2real_rec(vars, degs, var_R * realPart - var_I * imagPart, var_R * imagPart + var_I * realPart, cur_ind, Degree(cur_deg.explvar-1, cur_deg.conjvar))
			elseif cur_deg.conjvar > 0
				cur_deg.explvar == 0 || warn("cur_deg.explvar should be 0 (and not $(cur_deg.explvar)), set to this value")
				return cplx2real_rec(vars, degs, var_R * realPart + var_I * imagPart, var_R * imagPart - var_I * realPart, cur_ind, Degree(0, cur_deg.conjvar-1))
			end
		elseif isbool(var)
			return cplx2real_rec(vars, degs, var*realPart, var*imagPart, cur_ind, Degree(0,0))
		else
			return cplx2real_rec(vars, degs, Exponent(SortedDict(var=>cur_deg))*realPart, Exponent(SortedDict(var=>cur_deg))*imagPart, cur_ind, Degree(0,0))
		end
	end
end

"""
	realPart, imagPart = cplx2real(pol::Polynomial)

	Convert a complex polynomial in complex variables into `realPart` and
	`imagPart` polynomials of twice as many variables, real and imag parts of
	`pol` variables.
"""
function cplx2real(pol::Polynomial)
	realPart = Polynomial()
	imagPart = Polynomial()

	for (expo, λ) in pol
		realexpo, imagexpo = cplx2real(expo)

		realPart += realexpo*real(λ) - imagexpo*imag(λ)
		imagPart += imagexpo*real(λ) + realexpo*imag(λ)
	end
	return (realPart, imagPart)
end


function cplx2real_add(pol::Polynomial)
  realPart = Polynomial()
  imagPart = Polynomial()

  for (expo, λ) in pol
    realexpo, imagexpo = cplx2real(expo)

    add!(realPart, realexpo*real(λ) - imagexpo*imag(λ))
    add!(imagPart, imagexpo*real(λ) + realexpo*imag(λ))
  end
  return (realPart, imagPart)
end

function cplx2real(pt_C::Point)
	pt = Point()
	for (var, val) in pt_C
		if var.kind <: Complex
			var_real, var_imag = varname_cplx2real(var.name)
			if real(val) != 0
				pt[Variable(var_real, Real)] = real(val)
			end
			if imag(val) != 0
				pt[Variable(var_imag, Real)] = imag(val)
			end
		else
			if val != 0
				pt[var] = real(val)
			end
		end
	end
	return pt
end


function real2cplx(pt::Point)
	ptC = Point()
	for (var, val) in pt
		if ismatch(r"_Re$", var.name)
			var_c = Variable(var.name[1:end-3], Complex)
			if !haskey(ptC, var_c)
				ptC[var_c] = val
			else
				ptC[var_c] += val
			end

		elseif ismatch(r"_Im$", var.name)
			var_c = Variable(var.name[1:end-3], Complex)
			if !haskey(ptC, var_c)
				ptC[var_c] = val*im
			else
				ptC[var_c] += val*im
			end

		else
			@assert !iscomplex(var)
			ptC[var] = val
		end
	end
	return ptC
end
