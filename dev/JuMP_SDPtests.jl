using JuMP, Mosek
import MathProgComplex

m = Model(solver=MosekSolver())

@variable(m, Xs1[1:3,1:3], SDP)
@variable(m, Xs2[1:1,1:1], SDP)
@variable(m, Xs3[1:1,1:1], SDP)

# moment constraint one
@constraint(m, Xs1[1,1] == 1)

@constraint(m, Xs2[1,1] + Xs1[2,2] + Xs1[3,3] == 4)
@constraint(m, Xs3[1,1] - Xs1[2,1] - Xs1[3,1] == 0)

# Find upper bound
@objective(m, Max, -Xs1[1,3])

println(m)
solve(m)

MathProgBase.writeproblem(m.internalModel, "JuMP_task.jtask")
mv("JuMP_task.jtask", "JuMP_task.json", remove_destination = true)

println("Maximum value is ", getobjective(m))

println("\nSolution is:")
@show getvalue(Xs1)
@show getvalue(Xs2)
@show getvalue(Xs3)

@show getdual(Xs1)
@show getdual(Xs2)
@show getdual(Xs3)

############ Plain Mosek:

printstream(msg::String) = print(msg)

barvardim = [3, 1, 1]

numcon = 3
bkc = [Mosek.Boundkey(2), Mosek.Boundkey(2), Mosek.Boundkey(2)]
blc = [1, 4, 0]
buc = [1, 4, 0]

barai    = [1, 2, 2, 2, 3, 3, 3]
baraj    = [1, 1, 1, 2, 1, 1, 3]
barak    = [1, 2, 3, 1, 2, 3, 1]
baral    = [1, 2, 3, 1, 1, 1, 1]
baraijkl = [1, 1, 1, 1, 1/2, 1/2, 1/2]

barcj    = [1]
barck    = [2]
barcl    = [1]
barcjkl  = [-0.5]

# Create a task object and attach log stream printer
maketask() do task

    putstreamfunc(task, MSK_STREAM_LOG, printstream)

    # Append SDP matrix variables and scalar variables.
    # The variables will initially be fixed at zero.
    appendbarvars(task,barvardim)


    # Append 'numcon' empty constraints.
    # The constraints will initially have no bounds.
    appendcons(task,numcon)
    putconboundslice(task,1,numcon+1, bkc,blc,buc)

    # putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE)
    putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE)

    # Set constraints SDP vars coeffs
    putbarablocktriplet(task, length(barai), barai, baraj, barak, baral, baraijkl)

    # Objective matrices and constant
    putbarcblocktriplet(task, length(barcj), barcj, barck, barcl, barcjkl)

    writedata(task, "Mosek.jtask")
    mv("Mosek.jtask", "Mosek.json", remove_destination = true)

    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)

    info("suc:")
    @show getsuc(task, MSK_SOL_ITR)
    info("slc:")
    @show getslc(task, MSK_SOL_ITR)
    info("sux:")
    @show getsux(task, MSK_SOL_ITR)
    info("slx:")
    @show getslx(task, MSK_SOL_ITR)

    # Get status information about the solution
    prosta = getprosta(task, MSK_SOL_ITR)
    solsta = getsolsta(task, MSK_SOL_ITR)

    info("Primal solution:")
    @show getbarxj(task, MSK_SOL_ITR, 1)
    @show getbarxj(task, MSK_SOL_ITR, 2)

    info("Dual solution:")
    @show getbarsj(task, MSK_SOL_ITR, 1)
    @show getbarsj(task, MSK_SOL_ITR, 2)

end

# ##########################################################################
# m = Model(solver=MosekSolver())

# @variable(m, Xs1[1:6,1:6], SDP)
# @variable(m, Xs2[1:1,1:1], SDP)
# @variable(m, Xs3[1:1,1:1], SDP)

# # moment constraint one
# @constraint(m, Xs1[1,1] == 1)

# # Ball constraint
# @constraint(m, Xs1[4,4] + Xs1[6,6] + Xs2[1,1] == 16)
# @constraint(m, Xs1[2,2] + Xs1[3,3] + Xs3[1,1] == 100)

# # Find upper bound
# @objective(m, Max, Xs1[1,2])

# println(m)
# solve(m)
# println("Maximum value is ", getobjective(m))

# println("\nSolution is:")
# @show getvalue(Xs1)
# @show getvalue(Xs2)

# @show getdual(Xs1)
# @show getdual(Xs2)