using OPFInstances, BenchmarkTools

println("Loading module MathPrgoComplex")
@time using MathProgComplex

function main()
    instancepath = getinstancepath("Matpower", "QCQP", "case300")
    println("Working on $(splitdir(instancepath))")

    pb_c, pt = import_from_dat(instancepath)

    println("import_from_dat call:")
    @btime import_from_dat(instancepath);


    println("\npb_cplx2real call:")
    @btime pb_cplx2real(pb_c);

    return
end

main()
