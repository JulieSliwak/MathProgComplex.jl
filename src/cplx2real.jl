export cplx2real, pb_cplx2real, real2cplx


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


# Conversion of all complex variables to real ones in the Poly structures

"""
	realPart, imagPart = cplx2real(expo::Exponent)

	Convert a complex Exponent in complex variables into `realPart` and
	`imagPart` polynomials of twice as many variables, real and imag parts of
	`expo` variables. Done recursively with the `cplx2real_rec` function.
"""
function cplx2real(expo::Exponent, λ::Number)
	a=real(λ)
	b=imag(λ)
	if expo.degree.explvar  == 1 && expo.degree.conjvar  == 1
		if length(expo.expo)==1
			realPart = Polynomial(); add!(realPart, 0)
			imagPart = Polynomial(); add!(imagPart, 0)

			z1 = first(expo.expo)
			# println("Square : ", z1[1])
			x1 = Variable(z1[1].name*"_Re", Real)
			y1 = Variable(z1[1].name*"_Im", Real)
			x1x1 = SortedDict{Variable, Degree}(x1=>Degree(2,0))
			y1y1 = SortedDict{Variable, Degree}(y1=>Degree(2,0))
			# (a+ib) * conj(z1) * z1 = ax1^2+ay1^2+ibx1^2+iby1^2
			if a != 0
				realPart = Polynomial(SortedDict{Exponent, Number}(Exponent(x1x1)=>a, Exponent(y1y1)=>a))
			end
			if b != 0
				imagPart = Polynomial(SortedDict{Exponent, Number}(Exponent(x1x1)=>b, Exponent(y1y1)=>b))
			end
			return realPart, imagPart
		else
			if first(expo.expo)[2].conjvar==1
				z1 = first(expo.expo)
				z2 = last(expo.expo)
			else
				z2 = first(expo.expo)
				z1 = last(expo.expo)
			end
			# println("Bilinear : ", z1, ", ", z2)
			x1 = Variable(z1[1].name*"_Re", Real)
			y1 = Variable(z1[1].name*"_Im", Real)
			x2 = Variable(z2[1].name*"_Re", Real)
			y2 = Variable(z2[1].name*"_Im", Real)
			# (a+ib)*conj(z1) * z2  = (a+ib)*(x1-iy1)*(x2+iy2)
			# = (a+ib)*(x1.x2 + y1.y2 + ix1.y2-ix2.y1)
			# = a.x1.x2 + a.y1.y2 + ia.x1.y2-ia.x2.y1+ib.x1.x2 + ib.y1.y2 - b.x1.y2+b.x2.y1
			# =
			#    a.x1.x2+a.y1.y2-b.x1.y2+b.x2.y1
			#    +
			#    i(a.x1.y2-a.x2.y1+b.x1.x2+b.y1.y2)
			x1x2 = SortedDict{Variable, Degree}(x1=>Degree(1,0), x2=>Degree(1,0))
			y1y2 = SortedDict{Variable, Degree}(y1=>Degree(1,0), y2=>Degree(1,0))
			x1y2 = SortedDict{Variable, Degree}(x1=>Degree(1,0), y2=>Degree(1,0))
			x2y1 = SortedDict{Variable, Degree}(x2=>Degree(1,0), y1=>Degree(1,0))

			realPart = Polynomial(SortedDict{Exponent, Number}(Exponent(x1x2)=>a, Exponent(y1y2)=>a, Exponent(x1y2)=>-b, Exponent(x2y1)=>b))
			imagPart = Polynomial(SortedDict{Exponent, Number}(Exponent(x1x2)=>b, Exponent(y1y2)=>b, Exponent(x1y2)=>a, Exponent(x2y1)=>-a))
			return realPart, imagPart
		end
	elseif expo.degree.explvar  == 0 && expo.degree.conjvar  == 0
		realPart = Polynomial(); add!(realPart, a)
		imagPart = Polynomial(); add!(imagPart, b)
		return realPart, imagPart
	else
		realPart = Polynomial(); add!(realPart, 1)
		imagPart = Polynomial(); add!(imagPart, 0)

		vars, inds = collect(keys(expo.expo)), collect(values(expo.expo))
		return cplx2real_rec(vars, inds, realPart, imagPart, length(expo)+1, Degree(0,0))
	end
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
				realPart_new = product(var_R, realPart)
				add!(realPart_new, product(-1, product(var_I, imagPart)))

				imagPart_new = product(var_R, imagPart)
				add!(imagPart_new,  product(var_I, realPart))

				return cplx2real_rec(vars, degs, realPart_new, imagPart_new, cur_ind, Degree(cur_deg.explvar-1, cur_deg.conjvar))
			elseif cur_deg.conjvar > 0
				cur_deg.explvar == 0 || warn(LOGGER, "cur_deg.explvar should be 0 (and not $(cur_deg.explvar)), set to this value")
				realPart_new = product(var_R, realPart)
				add!(realPart_new,  product(var_I, imagPart))

				imagPart_new = product(var_R, imagPart)
				add!(imagPart_new, product(-1, product(var_I, realPart)))

				return cplx2real_rec(vars, degs, realPart_new, imagPart_new, cur_ind, Degree(0, cur_deg.conjvar-1))
			end
		elseif isbool(var)
			return cplx2real_rec(vars, degs, var*realPart, var*imagPart, cur_ind, Degree(0,0))
		else
			return cplx2real_rec(vars, degs, product(Exponent(SortedDict(var=>cur_deg)), realPart),
											 product(Exponent(SortedDict(var=>cur_deg)), imagPart),
											 cur_ind,
											 Degree(0,0))
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
	if (expo.degree.explvar  == 1 && expo.degree.conjvar  == 1 )||( expo.degree.explvar  == 0 && expo.degree.conjvar  == 0)
		realexpo, imagexpo = cplx2real(expo, λ)

		add!(realPart, realexpo)
		add!(imagPart, imagexpo)
	else
	    realexpo, imagexpo = cplx2real(expo, λ)

		add!(realPart, product(realexpo,   real(λ)))
		add!(realPart, product(imagexpo, - imag(λ)))
		add!(imagPart, product(imagexpo,   real(λ)))
		add!(imagPart, product(realexpo,   imag(λ)))
	end
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
		if occursin(r"_Re$", var.name)
			var_c = Variable(var.name[1:end-3], Complex)
			if !haskey(ptC, var_c)
				ptC[var_c] = val
			else
				ptC[var_c] += val
			end

		elseif occursin(r"_Im$", var.name)
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
