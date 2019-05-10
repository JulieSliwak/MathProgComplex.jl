using MathProgComplex, DataStructures, OPFInstances
using JuMP, KNITRO

function main()

    problem_c, pt = import_from_dat(getinstancepath("Matpower", "QCQP", "case9"))
    problem = pb_cplx2real(problem_c)

    solver = KnitroSolver(KTR_PARAM_OUTLEV=3,
                        KTR_PARAM_MAXIT=600,
                        KTR_PARAM_SCALE=0,
                        KTR_PARAM_FEASTOL=1.0,
                        KTR_PARAM_OPTTOL=1.0,
                        KTR_PARAM_FEASTOLABS=1e-6,
                        KTR_PARAM_OPTTOLABS=1e-3,
                        KTR_PARAM_BAR_INITPT=2,
                        KTR_PARAM_PRESOLVE=0,
                        KTR_PARAM_HONORBNDS=0)

    m, JuMPvars = get_JuMP_cartesian_model(problem, solver)
    solve(m)

    println("Objective value : $(getobjectivevalue(m))\n")

    for (varname, var) in JuMPvars
      value = getvalue(var)
      println("$varname   $value")
    end


    # amplexportpath = joinpath(Pkg.dir("MathProgComplex"), "Knitro_runs")
    # point = Point()

    # export_to_dat(problem, amplexportpath)
    # run_knitro(amplexportpath)

end

main()
