using MathProgComplex, DataStructures, OPFInstances, BenchmarkTools
using JuMP, KNITRO
# using ProfileView

function main()

    instancepath = getinstancepath("Matpower", "QCQP", "case1354pegase")

    println("import_from_dat time:")
    @btime import_from_dat($instancepath)
    problem_c, pt = import_from_dat(instancepath)

    println("pb_cplx2real time:")
    @btime problem = pb_cplx2real($problem_c)
    problem = pb_cplx2real(problem_c)


    @profile pb_cplx2real(problem_c)

    open("prof.txt", "w") do s
        Profile.print(IOContext(s, :displaysize => (24, 500)), format=:flat)
        # Profile.print(IOContext(s))
     end

    ProfileView.view()

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

    println("get_JuMP_cartesian_model time:")
    # @btime get_JuMP_cartesian_model($problem, $solver)
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

# Timings, case9
# Without Manuel's update to order.jl, Exponent.jl

# import_from_dat time:
#   9.342 ms (30396 allocations: 1.27 MiB)
# pb_cplx2real time:
#   38.879 ms (370727 allocations: 16.47 MiB)
# get_JuMP_cartesian_model time:
#   892.305 μs (8573 allocations: 429.27 KiB)

# import_from_dat time:
#   7.991 ms (18790 allocations: 756.09 KiB)
# pb_cplx2real time:
#   16.707 ms (136490 allocations: 5.73 MiB)
# get_JuMP_cartesian_model time:
#   874.254 μs (8574 allocations: 429.52 KiB)



#   case300

#   import_from_dat time:
#   331.475 ms (764052 allocations: 29.46 MiB)
# pb_cplx2real time:
#   930.399 ms (6244658 allocations: 255.18 MiB)
# get_JuMP_cartesian_model time:
#   52.450 ms (395145 allocations: 19.43 MiB)

#   import_from_dat time:
#   457.236 ms (1841742 allocations: 79.81 MiB)
# pb_cplx2real time:
#   2.460 s (20014923 allocations: 889.81 MiB)
# get_JuMP_cartesian_model time:
#   54.524 ms (395145 allocations: 19.43 MiB)

