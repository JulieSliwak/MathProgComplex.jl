using JuMP, KNITRO, MathProgComplex

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end


@testset "Matpower case9 KNITRO" begin
    instancepath = joinpath(Pkg.dir("MathProgComplex"), "test", "instances")
    pb_real, ~ = import_from_dat(joinpath(instancepath,"case9real.dat"))
    mysolver = KnitroSolver(KTR_PARAM_OUTLEV=3,
                          KTR_PARAM_MAXIT=600,
                          KTR_PARAM_SCALE=0,
                          KTR_PARAM_FEASTOL=1.0,
                          KTR_PARAM_OPTTOL=1.0,
                          KTR_PARAM_FEASTOLABS=1e-6,
                          KTR_PARAM_OPTTOLABS=1e-3,
                          KTR_PARAM_BAR_INITPT=2,
                          KTR_PARAM_PRESOLVE=0,
                          KTR_PARAM_HONORBNDS=0)

    m, variables_jump = get_JuMP_cartesian_model(pb_real, mysolver)
    solve(m)
    @test getobjectivevalue(m) ≈ 1.45883471040144e+003 atol=1e-6
end
