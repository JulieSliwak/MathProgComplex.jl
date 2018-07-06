using OPFInstances, BenchmarkTools, DataStructures

println("Loading module MathProgComplex")
tic();
using MathProgComplex
toc();

function main()
    instancepath = getinstancepath("Matpower", "QCQP", "case300")
    println("\nWorking on $(splitdir(instancepath))")

    seek_efficiency!(true)
    @show seek_efficiency()

    pb_c, pt = import_from_dat(instancepath)

    println("import_from_dat call:")
    @btime import_from_dat($instancepath);
    @show Base.summarysize(pb_c)

    pb = pb_cplx2real(pb_c)

    println("\npb_cplx2real call:")
    @btime pb_cplx2real($pb_c);
    @show Base.summarysize(pb)

    return
end

function main_work()
    instancepath = getinstancepath("Matpower", "QCQP", "case300")
    println("\nWorking on $(splitdir(instancepath))")

    seek_efficiency!(true)
    @show seek_efficiency()
    pb_c, pt = import_from_dat(instancepath)

    pb_cplx2real(pb_c)
    seek_efficiency!(false)

    return
end


# main_work()
main()
