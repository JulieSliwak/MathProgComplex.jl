using MathProgComplex


function build_WB2real()
    instancepath = joinpath(Pkg.dir("MathProgComplex"), "test", "instances")
    WB5_cplx, initpt = import_from_dat(instancepath, filename="WB5.dat")

    WB5 = pb_cplx2real(WB5_cplx)
    print(WB5)

    exportpath = joinpath(Pkg.dir("MathProgComplex"), "Knitro_runs")
    !ispath(exportpath) && mkpath(exportpath)

    export_to_dat(WB5, exportpath)

    run_knitro(exportpath)

    return 
end


build_WB2real()