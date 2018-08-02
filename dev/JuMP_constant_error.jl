using JuMP, Mosek, SCS
import MathProgComplex

function run_SDPpb(cst, optsense)
    m = Model(solver=MosekSolver(LOG=0))
    # m = Model(solver=SCSSolver(verbose=0))

    @variable(m, Xs1[1:3,1:3], SDP)
    @variable(m, Xs2[1:1,1:1], SDP)
    @variable(m, Xs3[1:1,1:1], SDP)

    @constraint(m, Xs1[1,1] == 1)
    @constraint(m, Xs2[1,1] + Xs1[2,2] + Xs1[3,3] == 4)
    @constraint(m, Xs3[1,1] - Xs1[2,1] - Xs1[3,1] == 0)

    # Find upper bound
    @objective(m, optsense, -Xs1[1,3]+cst)
    solve(m)
    return getobjectivevalue(m)
end

function main()
    @show run_SDPpb(0., :Min)
    @show run_SDPpb(0.5, :Min)

    @show run_SDPpb(0., :Max)
    @show run_SDPpb(0.5, :Max)

    return nothing
end

main()
