using MathProgComplex


function build_WB2real()
    WB2 = Problem()

    VOLT1_Im = Variable("VOLT1_Im", Real)
    VOLT1_Re = Variable("VOLT1_Re", Real)
    VOLT2_Im = Variable("VOLT2_Im", Real)
    VOLT2_Re = Variable("VOLT2_Re", Real)

    set_objective!(WB2, (192.3076923076923)*VOLT1_Im^2 + (-192.3076923076923)*VOLT1_Im * VOLT2_Im + (961.5384615384615)*VOLT1_Im * VOLT2_Re + (192.3076923076923)*VOLT1_Re^2 + (-961.5384615384615)*VOLT1_Re * VOLT2_Im + (-192.3076923076923)*VOLT1_Re * VOLT2_Re)

    add_constraint!(WB2, "BaseCase_1_BALANCE-UNIT_Im", -400.0 << ((480.7692307692308)*VOLT1_Im^2 + (-480.7692307692308)*VOLT1_Im * VOLT2_Im + (-96.15384615384615)*VOLT1_Im * VOLT2_Re + (480.7692307692308)*VOLT1_Re^2 + (96.15384615384615)*VOLT1_Re * VOLT2_Im + (-480.7692307692308)*VOLT1_Re * VOLT2_Re) << 400.0)
    add_constraint!(WB2, "BaseCase_1_BALANCE-UNIT_Re", 0.0 << ((96.15384615384615)*VOLT1_Im^2 + (-96.15384615384615)*VOLT1_Im * VOLT2_Im + (480.7692307692308)*VOLT1_Im * VOLT2_Re + (96.15384615384615)*VOLT1_Re^2 + (-480.7692307692308)*VOLT1_Re * VOLT2_Im + (-96.15384615384615)*VOLT1_Re * VOLT2_Re) << 600.0)
    add_constraint!(WB2, "BaseCase_1_Volt_VOLTM_Re", 0.9025 << (VOLT1_Im^2 + VOLT1_Re^2) << 1.1025)
    add_constraint!(WB2, "BaseCase_2_BALANCE-LOAD_Im", -350.0 + (-480.7692307692308)*VOLT1_Im * VOLT2_Im + (96.15384615384615)*VOLT1_Im * VOLT2_Re + (-96.15384615384615)*VOLT1_Re * VOLT2_Im + (-480.7692307692308)*VOLT1_Re * VOLT2_Re + (480.7692307692308)*VOLT2_Im^2 + (480.7692307692308)*VOLT2_Re^2 == 0.0)
    add_constraint!(WB2, "BaseCase_2_BALANCE-LOAD_Re", 350.0 + (-96.15384615384615)*VOLT1_Im * VOLT2_Im + (-480.7692307692308)*VOLT1_Im * VOLT2_Re + (480.7692307692308)*VOLT1_Re * VOLT2_Im + (-96.15384615384615)*VOLT1_Re * VOLT2_Re + (96.15384615384615)*VOLT2_Im^2 + (96.15384615384615)*VOLT2_Re^2 == 0.0)
    add_constraint!(WB2, "BaseCase_2_Volt_VOLTM_Re", 0.9025 << (VOLT2_Im^2 + VOLT2_Re^2) << 1.056784)

    exportpath = joinpath(Pkg.dir("MathProgComplex"), "Knitro_runs")
    !ispath(exportpath) && mkpath(exportpath)

    export_to_dat(WB2, exportpath)

    run_knitro(exportpath)

    return WB2
end


build_WB2real()