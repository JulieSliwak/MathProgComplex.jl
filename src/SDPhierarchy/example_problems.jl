export buildPOP_1v1c, buildPOP_1v2c, buildPOP_1v2, buildPOP_EllJoszMolc, buildPOP_knapsack
export lasserre_ex1, lasserre_ex2, lasserre_ex3, lasserre_ex5
export test_gloptipoly_ex1, test_gloptipoly_ex2
export get_WB5cliques, get_case9cliques

############################################
### Small toy geometric problems
############################################

function buildPOP_1v1c()
    z = Variable("z", Complex)
    problem = Problem()
    set_objective!(problem, -real(z))
    add_constraint!(problem, "ineq", abs2(z) << 4)
    return problem
end

function buildPOP_1v2c()
    z = Variable("z", Complex)
    problem = Problem()
    set_objective!(problem, imag(z))
    add_constraint!(problem, "ineq_brn", abs2(z) << 1)
    θ = π/3
    add_constraint!(problem, "ineq_rot", real(z*exp(-im*θ)) >> 0)
    return problem
end

function buildPOP_1v2()
    x1 = Variable("x", Real)
    x2 = Variable("y", Real)
    problem = Problem()
    set_objective!(problem, -1.0*x1)
    add_constraint!(problem, "ineq", (x1^2+x2^2) << 4)
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
    set_objective!(problem, 3-abs2(z1)-0.5*im*z1*conj(z2)^2+0.5im*z2^2*conj(z1))
    add_constraint!(problem, "eq1", (abs2(z1)-0.25*z1^2-0.25*conj(z1)^2) == 1)
    add_constraint!(problem, "eq2", (abs2(z1)+abs2(z2)) == 3)
    add_constraint!(problem, "eq3", (im*z2-im*conj(z2)) == 0)
    add_constraint!(problem, "ineq", (z2+conj(z2)) >> 0)
    return problem
end

function buildPOP_knapsack()
    x1 = Variable("x1", Bool)
    x2 = Variable("x2", Bool)
    x3 = Variable("x3", Bool)
    x4 = Variable("x4", Bool)
    x5 = Variable("x5", Bool)
    problem = Problem()
    set_objective!(problem, -1*(x1+x2+x3+x4+1.5*x5))         # Problems are minimization only right now
    add_constraint!(problem, "knapsack", (47*x1+45*x2+79*x3+53*x4+53*x5) << 178)
    return problem
end


############################################
### Global Optim pbs from Lasserre2001
############################################

"""
    problem, relax_ctx = lasserre_ex1()

From Lasserre2001, global minimum : (3) -0.2428.
"""
function lasserre_ex1()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    set_objective!(problem, (x1^2+1)^2 + (x2^2+1)^2 + (x1+x2+1)^2)

    ε_abs=1e-2
    order_to_obj = OrderedDict(2=>2.5074)

    return problem, order_to_obj, ε_abs
end

"""
    problem, relax_ctx = lasserre_ex2()

From Lasserre2001, global minimum : -11.4581.
"""
function lasserre_ex2()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    set_objective!(problem, (x1^2+1)^2 + (x2^2+1)^2 -2*(x1+x2+1)^2)

    ε_abs=1e-3
    order_to_obj = OrderedDict(3=>-11.4581)

    return problem, order_to_obj, ε_abs
end

"""
    problem, relax_ctx = lasserre_ex3()

From Lasserre2001, global minimum : -1/27, x1*² = x2*² = 1/3.
"""
function lasserre_ex3()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    set_objective!(problem, x1^2 * x2^2 * (x1^2 + x2^2 - 1))
    add_constraint(problem, "ball_constraint", (x1^2 + x2^2) << 10^2)

    ε_abs=1e-3
    order_to_obj = OrderedDict(3=>-1/27)

    return problem, order_to_obj, ε_abs
end

"""
    problem, relax_ctx = lasserre_ex5()

From Lasserre2001, global minimum : -2, for (1, 2).
Relaxation : order 1 -> -3; order 2 -> -2.
http://www.ii.uib.no/~lennart/drgrad/Lasserre2001.pdf
"""
function lasserre_ex5(;d = 2)
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    problem = Problem()
    set_objective!(problem, -(x1-1)^2 -(x1-x2)^2 -(x2-3)^2)
    add_constraint!(problem, "crt1", (1-(x1-1)^2) >> 0)
    add_constraint!(problem, "crt2", (1-(x1-x2)^2) >> 0)
    add_constraint!(problem, "crt3", (1-(x2-3)^2) >> 0)

    ε_abs=1e-3
    order_to_obj = OrderedDict(1=>-3.0000,
                            2=>-2.0000)

    return problem, order_to_obj, ε_abs
end


"""
    problem, order_to_obj, ε_abs = test_gloptipoly_ex1()

Return the problem, order to solution value and absolute tolerance
for the SDP test case taken from the gloptimoly doc:
http://homepages.laas.fr/henrion/papers/gpcocos.pdf p7

Order : 
  - 1 -> -6.0000
  - 2 -> -5.6923
  - 3 -> -4.0685
  - 3 -> -4.0000 (exact)
"""
function test_gloptipoly_ex1()
    problem = Problem()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    x3 = Variable("x3", Real)

    set_objective!(problem, -2*x1+x2-x3)
    add_constraint!(problem, "ctr1", (x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24) >> 0)
    add_constraint!(problem, "ctr2", (x1+x2+x3) << 4)
    add_constraint!(problem, "ctr3", (3*x2+x3) << 6)
    add_constraint!(problem, "def_x1", 0 << x1 << 2)
    add_constraint!(problem, "def_x2", 0 << x2)
    add_constraint!(problem, "def_x3", 0 << x3 << 3)

    ε_abs=1e-3
    order_to_obj = OrderedDict(1=>-6.0000,
                            2=>-5.6923,
                            3=>-4.0685)

    return problem, order_to_obj, ε_abs
end

"""
    problem, order_to_obj, ε_abs = test_gloptipoly_ex2()

Return the problem, order to solution value and absolute tolerance
for the SDP test case taken from the gloptimoly doc:
http://homepages.laas.fr/henrion/papers/gloptipoly.pdf p14

Order : 
  - 1 -> primal infeasible
  - 2 -> 17.9189
  - 3 -> 17.0000 (exact)
"""
function test_gloptipoly_ex2()
    problem = Problem()
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    x3 = Variable("x3", Real)
    x4 = Variable("x4", Real)
    x5 = Variable("x5", Real)

    set_objective!(problem, 42*x1+44*x2+45*x3+47*x4+47.5*x5-50*(x1^2+x2^2+x3^2+x4^2+x5^2))
    add_constraint!(problem, "ctr1", (20*x1+12*x2+11*x3+7*x4+4*x5) << 40)
    add_constraint!(problem, "def_x1", 0 << x1 << 1)
    add_constraint!(problem, "def_x2", 0 << x2 << 1)
    add_constraint!(problem, "def_x3", 0 << x3 << 1)
    add_constraint!(problem, "def_x4", 0 << x4 << 1)
    add_constraint!(problem, "def_x5", 0 << x5 << 1)

    ε_abs=1e-3
    order_to_obj = SortedDict(2=>-17.9189)      # order 1 relaxation is primal infeasable

    return problem, order_to_obj, ε_abs
end



############################################
### Clique decomposition of small OPF instances
############################################

function get_WB5cliques(relax_ctx, problem)
    if !relax_ctx.issparse
        return get_maxcliques(relax_ctx, problem)
    else
        maxcliques = DictType{String, Set{Variable}}()
        maxcliques["clique1"] = Set{Variable}([
            Variable("VOLT_1_Im", Real),
            Variable("VOLT_1_Re", Real),
            Variable("VOLT_2_Im", Real),
            Variable("VOLT_2_Re", Real),
            Variable("VOLT_3_Im", Real),
            Variable("VOLT_3_Re", Real)])
        maxcliques["clique2"] = Set{Variable}([
            Variable("VOLT_2_Im", Real),
            Variable("VOLT_2_Re", Real),
            Variable("VOLT_3_Im", Real),
            Variable("VOLT_3_Re", Real),
            Variable("VOLT_4_Im", Real),
            Variable("VOLT_4_Re", Real),
            Variable("VOLT_5_Im", Real),
            Variable("VOLT_5_Re", Real)])
        return maxcliques
    end
end

function get_case9cliques(relax_ctx, problem)
    if !relax_ctx.issparse
        return get_maxcliques(relax_ctx, problem)
    else
        maxcliques = DictType{String, Set{Variable}}()
        maxcliques["clique1"] = Set{Variable}([
            Variable("VOLT_1_Im", Real),
            Variable("VOLT_1_Re", Real),
            Variable("VOLT_5_Im", Real),
            Variable("VOLT_5_Re", Real),
            Variable("VOLT_4_Im", Real),
            Variable("VOLT_4_Re", Real),
            Variable("VOLT_9_Im", Real),
            Variable("VOLT_9_Re", Real),
            Variable("VOLT_8_Im", Real),
            Variable("VOLT_8_Re", Real)])
        maxcliques["clique2"] = Set{Variable}([
            Variable("VOLT_2_Im", Real),
            Variable("VOLT_2_Re", Real),
            Variable("VOLT_9_Im", Real),
            Variable("VOLT_9_Re", Real),
            Variable("VOLT_8_Im", Real),
            Variable("VOLT_8_Re", Real),
            Variable("VOLT_7_Im", Real),
            Variable("VOLT_7_Re", Real),
            Variable("VOLT_6_Im", Real),
            Variable("VOLT_6_Re", Real)])
        maxcliques["clique3"] = Set{Variable}([
            Variable("VOLT_3_Im", Real),
            Variable("VOLT_3_Re", Real),
            Variable("VOLT_7_Im", Real),
            Variable("VOLT_7_Re", Real),
            Variable("VOLT_6_Im", Real),
            Variable("VOLT_6_Re", Real),
            Variable("VOLT_5_Im", Real),
            Variable("VOLT_5_Re", Real),
            Variable("VOLT_4_Im", Real),
            Variable("VOLT_4_Re", Real)])
        return maxcliques
    end
end
