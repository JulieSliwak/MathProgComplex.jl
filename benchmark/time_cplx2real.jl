using OPFInstances, BenchmarkTools, DataStructures, Memento

println("Loading module MathProgComplex")
tic();
using MathProgComplex
toc();

function main(instance = "case300")
    instancepath = getinstancepath("Matpower", "QCQP", instance)
    println("\nWorking on $(splitdir(instancepath))")

    setlevel!(getlogger(MathProgComplex), "debug")

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

    setlevel!(getlogger(MathProgComplex), "debug")
    pb_c, pt = import_from_dat(instancepath)

    pb_cplx2real(pb_c)

    # turn off debug info
    setlevel!(getlogger(MathProgComplex), "info")

    return
end


# main_work()
# main("case300")
# main("case1354pegase")
