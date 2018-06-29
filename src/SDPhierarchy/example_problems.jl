include(joinpath(ROOT, "src_PowSysMod", "PowSysMod_body.jl"))


############################
### Geometric problems
############################

function buildPOP_1v1c()
    z = Variable("z", Complex)
    problem = Problem()
    add_variable!(problem, z)
    set_objective!(problem, -real(z))
    add_constraint!(problem, "ineq", abs2(z) << 4)
    return problem
end

function buildPOP_1v2c()
    z = Variable("z", Complex)
    problem = Problem()
    add_variable!(problem, z)
    set_objective!(problem, imag(z))
    add_constraint!(problem, "ineq_brn", abs2(z) << 1)
    θ = π/3
    add_constraint!(problem, "ineq_rot", real(z*exp(-im*θ)) >> 0)
    return problem
end

function buildPOP_1v2()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    add_variable!(problem, x1); add_variable!(problem, x2)
    set_objective!(problem, -1.0*x1)
    add_constraint!(problem, "ineq", (x1^2+x2^2) << 1)
    θ1 = π/3
    add_constraint!(problem, "eq_rot1", (cos(θ1)*x1+sin(θ1)*x2) == 0)
    # θ2 = -π/3
    # add_constraint!(problem, "ineq_rot2", (cos(θ2)*x1+sin(θ2)*x2) >> 0)
    return problem
end

"""
    problem = buildPOP_2v3c

    Elliptic example problemp from Josz, Molzahn 2018 paper.
"""
function buildPOP_EllJoszMolc()
    z1 = Variable("z1", Complex)
    z2 = Variable("z2", Complex)
    problem = Problem()
    add_variable!(problem, z1); add_variable!(problem, z2);
    set_objective!(problem, 3-abs2(z1)-0.5*im*z1*conj(z2)^2+0.5im*z2^2*conj(z1))
    add_constraint!(problem, "eq1", (abs2(z1)-0.25*z1^2-0.25*conj(z1)^2) == 1)
    add_constraint!(problem, "eq2", (abs2(z1)+abs2(z2)) == 3)
    add_constraint!(problem, "eq3", (im*z2-im*conj(z2)) == 0)
    add_constraint!(problem, "ineq", (z2+conj(z2)) >> 0)
    return problem
end

############################
### OPF problems
############################

function buildPOP_WB2(; v2max = 0.976, rmeqs = false, setnetworkphase=false, addball=false)
    OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", "WB2.m"))
    problem_c = build_globalpb!(OPFpbs)

    ## Converting to real ineq. only problem
    !rmeqs || change_eq_to_ineq!(problem_c)
    problem = pb_cplx2real(problem_c)

    if setnetworkphase
        ## Fixing volt phase of last bus to 0
        lastctr = problem.constraints["BaseCase_2_Volt_VOLTM_Re"]
        rm_constraint!(problem, "BaseCase_2_Volt_VOLTM_Re")

        ## Setting imag part to 0
        # pt = Point(SortedDict(Variable("BaseCase_2_VOLT_Im", Real)=>0.0), isdense=true)
        # infer_problem!(problem, pt)

        add_constraint!(problem, "BaseCase_2_Volt_VOLTM_Re", sqrt(lastctr.lb) << Variable("BaseCase_2_VOLT_Re", Real) << v2max)
        add_constraint!(problem, "BaseCase_2_Volt_VOLTM_Im", Variable("BaseCase_2_VOLT_Im", Real) == 0)
    else
        problem.constraints["BaseCase_2_Volt_VOLTM_Re"].ub = v2max^2
    end

    ## Adding ball constraint
    if addball
        p = Polynomial()
        for var in problem.variables
            p += Variable(var[1], var[2])^2
        end
        ub = problem.constraints["BaseCase_1_Volt_VOLTM_Re"].ub + v2max^2
        add_constraint!(problem, "Ball_ctr", p << ub)
    end

    return problem
end

function buildPOP_WB5(; q5min = 1.05, rmeqs = false)
    OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", "WB5.m"))
    Sgen = OPFpbs["BaseCase"].ds.bus["BUS_5"]["Gen_1"].power_min
    OPFpbs["BaseCase"].ds.bus["BUS_5"]["Gen_1"].power_min = real(Sgen) + im*q5min
    problem_c = build_globalpb!(OPFpbs)

    ## Converting to real ineq. only problem
    !rmeqs || change_eq_to_ineq!(problem_c)
    return pb_cplx2real(problem_c)
end



############################
### Global Optim pbs from Lasserre2001
############################

"""
    problem, relax_ctx = lasserre_ex1()

    From Lasserre2001, global minimum : (3) -0.2428.
"""
function lasserre_ex1()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    add_variable!(problem, x1); add_variable!(problem, x2)
    set_objective!(problem, (x1^2+1)^2 + (x2^2+1)^2 + (x1+x2+1)^2)

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = 2)
    return problem, relax_ctx
end

"""
    problem, relax_ctx = lasserre_ex2()

    From Lasserre2001, global minimum : -11.4581.
"""
function lasserre_ex2()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    add_variable!(problem, x1); add_variable!(problem, x2)
    set_objective!(problem, (x1^2+1)^2 + (x2^2+1)^2 -2*(x1+x2+1)^2)

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = 2)
    return problem, relax_ctx
end

"""
    problem, relax_ctx = lasserre_ex3()

    From Lasserre2001, global minimum : -1/27, x1*² = x2*² = 1/3.
"""
function lasserre_ex3()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    add_variable!(problem, x1); add_variable!(problem, x2)
    set_objective!(problem, x1^2 * x2^2 * (x1^2 + x2^2 - 1))

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = 3)
    return problem, relax_ctx
end

"""
    problem, relax_ctx = lasserre_ex5()

    From Lasserre2001, global minimum : -2, for (1, 2).
    Relaxation : order 1 -> -3; order 2 -> -2.
"""
function lasserre_ex5(;d = 2)
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    add_variable!(problem, x1); add_variable!(problem, x2)
    set_objective!(problem, -(x1-1)^2 -(x1-x2)^2 -(x2-3)^2)
    add_constraint!(problem, "crt1", (1-(x1-1)^2) >> 0)
    add_constraint!(problem, "crt2", (1-(x1-x2)^2) >> 0)
    add_constraint!(problem, "crt3", (1-(x2-3)^2) >> 0)

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = d)
    return problem, relax_ctx
end
