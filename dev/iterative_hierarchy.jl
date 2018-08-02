using DataStructures, SCS, Mosek, OPFInstances
using MathProgComplex

!isdefined(:MPC) && (const MPC = MathProgComplex)

function main()
    problem_c, pt = MPC.import_from_dat(getinstancepath("Matpower", "QCQP", "WB5"))
    problem = MPC.pb_cplx2real(problem_c)

    workpath = joinpath("Mosek_runs", "worksdp")
    ispath(workpath) && rm(workpath, recursive=true)
    mkpath(workpath)

    # cstobj = problem.objective[MPC.Exponent()]
    # cstobj = 0
    # problem.objective += -cstobj
    # @assert !haskey(problem.objective, MPC.Exponent())

    di = Dict("LOAD_2_Im" => 1,
            "LOAD_2_Re" => 1,
            "LOAD_3_Im" => 1,
            "LOAD_3_Re" => 1,
            "LOAD_4_Im" => 1,
            "LOAD_4_Re" => 1,
            "UNIT_1_Im" => 1,
            "UNIT_1_Re" => 1,
            "UNIT_5_Im" => 1,
            "UNIT_5_Re" => 1,
            "VOLTM_1_Re" => 1,
            "VOLTM_2_Re" => 1,
            "VOLTM_3_Re" => 1,
            "VOLTM_4_Re" => 1,
            "VOLTM_5_Re" => 1)

    relax_ctx = MPC.set_relaxation(problem; hierarchykind=:Real,
                                            issparse = false,
                                            d = 2,
                                            # di = di,
                                            params = Dict(:opt_outlev=>2,
                                                          :opt_outmode=>0,
                                                          :opt_pbsolved=>:MomentRelaxation,
                                                          :opt_solver=>:MosekCAPI))

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj_JuMP, dualobj_JuMP = MPC.run_hierarchy(problem, relax_ctx, indentedprint=true, save_pbs=true,
                                                                        primsol = primal,
                                                                        dualsol = dual)

    point1 = Point()
    for (varname, varkind) in problem.variables
        setindex!(point1, primal["clique1", varname, "1"], Variable(varname, varkind))
    end

    @show point1

    point2 = Point()
    for (varname, varkind) in problem.variables
        setindex!(point2, sqrt(primal["clique1", varname, varname]), Variable(varname, varkind))
    end

    @show point2

    println("slacks point:\n", get_slacks(problem, point1))

    println("slacks point2:\n", get_slacks(problem, point2))

    return nothing
end

main()
