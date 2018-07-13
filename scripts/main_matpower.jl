
using MathProgComplex, DataStructures, OPFInstances
using JuMP, KNITRO


function main()
    tic()
    problem_c, pt = import_from_dat(getinstancepath("Matpower", "QCQP", "case1354pegase"))
    #problem_c, pt = import_from_dat(getinstancepath("Matpower", "QCQP", "case9"))
    println("import_from_dat : ", toq())
    tic()
    @profile pb_cplx2real(problem_c)
    open("prof.txt", "w") do s
         Profile.print(IOContext(s, :displaysize => (24, 500)), format=:flat)
         # Profile.print(IOContext(s))
     end
    println("pb_cplx2real : ", toq())
    # tic()
    #
    # solver = KnitroSolver(KTR_PARAM_OUTLEV=3,
    #                     KTR_PARAM_MAXIT=600,
    #                     KTR_PARAM_SCALE=0,
    #                     KTR_PARAM_FEASTOL=1.0,
    #                     KTR_PARAM_OPTTOL=1.0,
    #                     KTR_PARAM_FEASTOLABS=1e-6,
    #                     KTR_PARAM_OPTTOLABS=1e-3,
    #                     KTR_PARAM_BAR_INITPT=2,
    #                     KTR_PARAM_PRESOLVE=0,
    #                     KTR_PARAM_HONORBNDS=0)
    #
    # println("KnitroSolver : ", toq())
    # tic()
    # m, JuMPvars = get_JuMP_cartesian_model(problem, solver)
    # println("get_JuMP_cartesian_model : ", toq())
    # tic()
    # solve(m)
    # println("solve : ", toq())
    #
    # println("Objective value : $(getobjectivevalue(m))\n")

   # for (varname, var) in JuMPvars
   #   value = getvalue(var)
   #   println("$varname   $value")
   # end
   #
   #
   #  amplexportpath = joinpath(Pkg.dir("MathProgComplex"), "Knitro_runs")
   #  point = Point()
   #
   #  export_to_dat(problem, amplexportpath)
   #  run_knitro(amplexportpath)

end

main()
