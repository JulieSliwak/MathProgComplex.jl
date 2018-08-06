export get_JuMP_cartesian_model, get_JuMP_solution


"""
    get_JuMP_cartesian_model(problem_poly::Problem, mysolver)

Return JuMP cartesian model associated to `mysolver` defined by `problem_poly`, a polynomial problem in complex numbers

# Example
```jldoctest
V1 = Variable("VOLT_1",Complex)
V2 = Variable("VOLT_2",Complex)
p_obj = V1*conj(V2) + (2+im)*abs2(V1) + 1+2im
p_ctr1 = abs2(V1)
p_ctr2 = 3im * V1*conj(V2) + abs2(V1)

problem_poly=Problem()
add_variable!(problem_poly,V1)
add_variable!(problem_poly,V2)
add_constraint!(problem_poly, "ctr1", 0.95^2 << p_ctr1 << 1.05^2)
add_constraint!(problem_poly, "ctr2", p_ctr2==0)
set_objective!(problem_poly, p_obj)
print(problem_poly)
mysolver = ClpSolver()
m = get_JuMP_cartesian_model(problem_poly, mysolver)
print(m)
```
"""
function get_JuMP_cartesian_model(problem_poly::Problem, mysolver)
    pb_poly_real = problem_poly
    m = Model(solver = mysolver)
    variables_jump = SortedDict{String, JuMP.Variable}()

    ## Define JuMP variables
    startReal = 1.1
    startBool = 0
    for (varname, vartype) in pb_poly_real.variables
    #     if init_point == Point()
    #         startReal = 1.1
    #         startBool = 0
    #     else
    #         startReal = init_point[Variable(varname, vartype)]
    #         startBool = init_point[Variable(varname, vartype)]
    #     end
        if vartype==Real
            variables_jump["$varname"] = @variable(m, basename="$varname", start= startReal)
        elseif vartype==Bool
            variables_jump["$varname"] = @variable(m, category=:Bin, basename="$varname", start=startBool)
        else
            error("$varname must be of type Real or Bool, here type is $vartype")
        end
    end

    ## Define JuMP constraints
    ctr_jump = SortedDict{String,JuMP.ConstraintRef}()
    ctr_exp = SortedDict{String,Tuple{Any,Float64,Float64}}()
    for (ctrname, modeler_ctr) in pb_poly_real.constraints
        polynome = modeler_ctr.p

        modeler_ctr.lb == -Inf-Inf*im || imag(modeler_ctr.lb) == 0 || error("get_JuMP_cartesian_model(): Found constraint with complex bound.\nOnly real models are valid at this point.")
        modeler_ctr.ub == +Inf+Inf*im || imag(modeler_ctr.ub) == 0 || error("get_JuMP_cartesian_model(): Found constraint with complex bound.\nOnly real models are valid at this point.")
        lb = real(modeler_ctr.lb)
        ub = real(modeler_ctr.ub)

        precond = modeler_ctr.precond
        for value in values(polynome)
            if imag(value)!=0
                error("Polynom coefficients have to be real numbers")
            end
        end
        my_timer = @elapsed s_ctr, ispolylinear = poly_to_NLexpression(m, variables_jump,polynome)
        ctr_exp[ctrname] = (s_ctr, lb, ub)
        # @printf("%-35s%10.6f s\n", "poly_to_NLexpression for $ctr", my_timer)
        if ispolylinear
            @constraint(m, lb <= s_ctr <= ub)
        else
            if precond == :sqrt
                if lb == -Inf && ub >0
                    ##TODO: sqrt constraints for Smax
                    ctr_jump[ctrname] = @NLconstraint(m, -Inf <= s_ctr <= ub)
                else
                    error("sqrt non applicable to $lb <= $s_ctr <= $ub")
                end
            else
                ctr_jump[ctrname] = @NLconstraint(m, lb <= s_ctr <= ub)
            end
        end
    end

    ## Define JuMP objective
    polynome_obj = pb_poly_real.objective
    my_timer = @elapsed s_obj, ispolylinear = poly_to_NLexpression(m, variables_jump,polynome_obj)
    # @printf("%-35s%10.6f s\n", "poly_to_NLexpression for objective", my_timer)
    if ispolylinear
        @objective(m, Min, s_obj)
    else
        @NLobjective(m,Min,s_obj)
    end
    return m, variables_jump, ctr_jump, ctr_exp
end

function poly_to_NLexpression(m::JuMP.Model, variables_jump::SortedDict{String, JuMP.Variable},polynome::Polynomial)
    s = 0
    ispolylinear = true
    d = compute_degree(polynome) #max degree on one variable, not total degree
    for expo in keys(polynome.poly)
        explsum, conjsum = get_sumdegs(expo)
        d.explvar = max(d.explvar, explsum)
        d.conjvar = max(d.conjvar,conjsum)
    end
    if d.explvar + d.conjvar > 1
        ispolylinear = false
    end
    for (monome,coeff) in polynome.poly
        prod = 1
        for (varname, degree) in monome.expo
            if !ispolylinear
                if degree.explvar > 1
                prod = @NLexpression(m, prod * variables_jump["$varname"]^degree.explvar)
                elseif degree.explvar == 1
                prod = @NLexpression(m, prod * variables_jump["$varname"])
                end
            else
                if degree.explvar != 1
                    error("Exponent is supposed to be degree 1 since polynome is linear. ")
                end
                prod = @expression(m, variables_jump["$varname"])
            end
        end
        if ispolylinear
            s = @expression(m, s + coeff * prod)
        else
            s = @NLexpression(m, s + coeff * prod)
        end
    end
    return s, ispolylinear
end


"""
    sol = get_JuMP_solution(m, variables_jump, rosenbrock)

    Return a `Point` storing the candidate point in `m`.
"""
function get_JuMP_solution(m, variables_jump, pb::Problem)
    sol = Point()
    for (varname, vartype) in pb.variables
        sol[MathProgComplex.Variable(varname, vartype)] = convert(vartype, getvalue(variables_jump[varname]))
    end
    return sol
end

###TEST
#
# ###########################################
# ###polynomial problem with complex variables
# V1 = Variable("VOLT_1",Complex)
# p_obj = 0.5*(V1+conj(V1))-0.5*im*(V1-conj(V1))
# p_ctr1 = abs2(V1)
# problem_poly=Problem()
# add_variable!(problem_poly,V1)
# add_constraint!(problem_poly, "ctr1", 1 << p_ctr1 << 1 )
# set_objective!(problem_poly, p_obj)
# print(problem_poly)
# ###########################################
# pb_poly_real = pb_cplx2real(problem_poly)
# println(pb_poly_real)
#
# mysolver = KnitroSolver(KTR_PARAM_MIP_INTVAR_STRATEGY=0, KTR_PARAM_OUTLEV=4, #=KTR_PARAM_HESSOPT=2=#)
# m = get_JuMP_cartesian_model(problem_poly, mysolver)
# print(m)
# solve(m)
#
#
# #################################################
# ## the same problem with JuMP
# mj = Model(solver=mysolver)
# @variable(mj, v1_re)
# @variable(mj, v1_im)
# @objective(mj, Min, v1_re+v1_im)
# @NLconstraint(mj, 1<= 1*v1_re^2 + 1*v1_im^2 <=1)
# print(mj)
# solve(mj)
# #################################################
#
# ###########################################
# ###polynomial problem with real variables
# V1_Re = Variable("VOLT_1_Re",Real)
# V1_Im = Variable("VOLT_1_Im",Real)
# p_obj = V1_Re + V1_Im
# p_ctr1 = V1_Re^2 + V1_Im^2
# problem_poly=Problem()
# add_variable!(problem_poly,V1_Re)
# add_variable!(problem_poly,V1_Im)
# add_constraint!(problem_poly, "ctr1", 1 << p_ctr1 << 1 )
# set_objective!(problem_poly, p_obj)
# print(problem_poly)
# mysolver = KnitroSolver(KTR_PARAM_MIP_INTVAR_STRATEGY=0, KTR_PARAM_OUTLEV=4)
# m = get_JuMP_cartesian_model(problem_poly, mysolver)
# print(m)
# solve(m)
#
# ###########################################
#






"""
    get_JuMP_polar_model(pb::Problem, mysolver)

Return JuMP polar model associated to `mysolver` defined by `pb`, a polynomial problem in complex numbers

# Example
```jldoctest
V1 = Variable("VOLT_1",Complex)
V2 = Variable("VOLT_2",Complex)
p_obj = V1*conj(V2) + (2+im)*abs2(V1) + 1+2im
p_ctr1 = abs2(V1)
p_ctr2 = Polynomial((3+4im) * V1*conj(V2))
pb=Problem()
add_variable!(pb,V1)
add_variable!(pb,V2)
add_constraint!(pb, "ctr1", 0 << p_ctr1 << 100)
# add_constraint!(pb, "ctr2", p_ctr2==3+4im)
set_objective!(pb, p_obj)
print(pb)
mysolver = IpoptSolver()
m = get_JuMP_polar_model(pb, mysolver)
print(m)
```
"""
function get_JuMP_polar_model(pb::Problem, mysolver)
    m = Model(solver = mysolver)
    jump_vars = SortedDict{String, JuMP.Variable}()
    for (varname, vartype) in pb.variables
        mod = "ρ_$varname"
        jump_vars["ρ_$varname"] = @variable(m, basename="ρ_$varname")
        jump_vars["θ_$varname"] = @variable(m, basename="θ_$varname")
    end
    ctr_jump = SortedDict{String,JuMP.ConstraintRef}()
    for (ctr, modeler_ctr) in pb.constraints
        p = modeler_ctr.p
        lb = modeler_ctr.lb
        ub = modeler_ctr.ub
        ps = SortedDict{Exponent, Tuple{Real, Real}}(exp=>(real(coeff), imag(coeff)) for (exp, coeff) in p)
        ## Real part of constraint
        @NLconstraint(m, real(lb) <= sum( coeff[1] * prod(jump_vars["ρ_$var"]^(exp[1]+exp[2]) for (var,exp) in mon) * cos(sum( (exp[1]-exp[2])*jump_vars["θ_$var"] for (var,exp) in mon))
                                        - coeff[2] * prod(jump_vars["ρ_$var"]^(exp[1]+exp[2]) for (var,exp) in mon) * sin(sum( (exp[1]-exp[2])*jump_vars["θ_$var"] for (var,exp) in mon))
                                        for (mon, coeff) in ps) <= real(ub))

        ## Imag part of constraint
        @NLconstraint(m, imag(lb) <= sum( coeff[2] * prod(jump_vars["ρ_$var"]^(exp[1]+exp[2]) for (var,exp) in mon) * cos(sum( (exp[1]-exp[2])*jump_vars["θ_$var"] for (var,exp) in mon))
                                        + coeff[1] * prod(jump_vars["ρ_$var"]^(exp[1]+exp[2]) for (var,exp) in mon) * sin(sum( (exp[1]-exp[2])*jump_vars["θ_$var"] for (var,exp) in mon))
                                        for (mon, coeff) in ps) <= imag(ub))

    polynome_obj = pb.objective
    ps = SortedDict{Exponent, Tuple{Real, Real}}(exp=>(real(coeff), imag(coeff)) for (exp, coeff) in polynome_obj)
    lb_real, ub_real = real(lb), real(ub)
    @NLobjective(m, :Min, lb_real <= sum( coeff[1] * prod(jump_vars["ρ_$var"]^(exp[1]+exp[2]) for (var,exp) in mon) * cos(sum( (exp[1]-exp[2])*jump_vars["θ_$var"] for (var,exp) in mon))
                                        - coeff[2] * prod(jump_vars["ρ_$var"]^(exp[1]+exp[2]) for (var,exp) in mon) * sin(sum( (exp[1]-exp[2])*jump_vars["θ_$var"] for (var,exp) in mon))
                                        for (mon, coeff) in ps) <= ub_real)
    end
    return m
end


####TEST
# V1 = Variable("VOLT_1",Complex)
# V2 = Variable("VOLT_2",Complex)
# p_obj = V1*conj(V2) + (2+im)*abs2(V1) + 1+2im
# p_ctr1 = abs2(V1)
# p_ctr2 = Polynomial((3+4im) * V1*conj(V2))
# pb=Problem()
# add_variable!(pb,V1)
# add_variable!(pb,V2)
# add_constraint!(pb, "ctr1", 0 << p_ctr1 << 100)
# # add_constraint!(pb, "ctr2", p_ctr2==3+4im)
# set_objective!(pb, p_obj)
# print(pb)
# mysolver = IpoptSolver()
# m = get_JuMP_polar_model(pb, mysolver)
# print(m)
