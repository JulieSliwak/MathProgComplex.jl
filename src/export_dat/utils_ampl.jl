# """
#   pt_sol = get_knitro_sol(filename)
#
# Compute a local optima of the filename complex QCQP.
# """
# function get_knitro_sol(filename::String)
#   open(joinpath("src_ampl", "minlp.dat"), "w") do f
#     println(f, "param: KEYS: LEFT RIGHT := include \"$(joinpath("..", "$filename"))\";")
#   end
#
#   cd("src_ampl")
#
#   instancename = split(filename, "\\")[end]
#   instancename = split(instancename, ".")[1]
#   out_log = "Knitro_$(instancename).log"
#   run(`cmd /c ampl minlp.run '>' $out_log`)
#   cd("..")
#
#   sol_cplx = read_Knitro_output("Knitro_sol.txt")
#   return sol_cplx
# end


"""
sol = read_Knitro_output(filepath)

Return the point corresponding to Knitro output stored in the filepath text
file.
"""
function read_Knitro_output(filepath::String)
  filepath = joinpath("src_ampl", "Knitro_sol.txt")
  data = readdlm(filepath)

  sol_cplx = Point()

  for i=1:size(data, 1)
    sol_cplx[Variable(data[i, 1], Complex)] = data[i, 2] + im*data[i, 3]
  end

  return sol_cplx
end



"""
  run_knitro(pb_path::String, src_ampl_path::String)

Run knitro on the `real_minlp.run`, `real_minlp.run` template scripts from the
`src_ampl_path` folder and `pb_path` files ("real_minlp_instance.dat" and
"real_minlp_precond_cstrs.dat").
"""
function run_knitro(pb_path::String, src_ampl_path::String)
    root = pwd()
    date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
    outlog = "Knitro_$(date).log"

    cd(pb_path)

    open("real_minlp.run", "w") do f
      println(f, "include $(joinpath(src_ampl_path, "real_minlp.run"));")
    end

    open("real_minlp.mod", "w") do f
      println(f, "include $(joinpath(src_ampl_path, "real_minlp.mod"));")
    end

    try
        run(`cmd /c ampl real_minlp.run '>' $(outlog)`)
        # run(`cmd /c ampl real_minlp.run `)
    catch
        warn("AMPL/Knitro failed, returning.")
    end

    cd(root)
end

"""
pt_knitro, pt_GOC = read_Knitro_output(pb_path, pb)

Read knitro output files at `pb_path` and return points according to variables
types in pb.
"""
function read_Knitro_output(pb_path::String, pb::Problem)
  function build_pt(vararray, valarray)
    pt = Point()
    for i=1:length(vararray)
      varname, vartype = vararray[i], get_variabletype(pb, String(vararray[i]))
      val = valarray[i]
      if (vartype == Bool)
        val = round(val)
      end
      pt[Variable(varname, vartype)] = val
    end
    return pt
  end

  root = pwd()
  cd(pb_path)

  ## Read Knitro solution
  files = filter(x->ismatch(r"solution\S+csv", x), readdir(pwd()))
  length(files) == 1 || warn("get_knitro_solutions(): $(length(files)) .csv files found in $(pb_path).\nExtracting solution from $(files[1]).")

  sol_data = readdlm(files[1], ';')

  pt_knitro = build_pt(sol_data[:, 1], sol_data[:, 2])
  pt_GOC = build_pt(sol_data[:, 1], sol_data[:, 3])

  cd(root)
  return pt_knitro, pt_GOC
end

function read_knitro_info_csvfile(filepath::String)
  lines = readdlm(joinpath(filepath, "knitro_info.csv"), ';')
  solve_result_1 = lines[2,2]
  opterror1 = lines[3,2]
  solve_result_2 = lines[5,2]
  opterror2 = lines[6,2]
  solve_result_3 = lines[8,2]
  opterror3 = lines[9,2]

  return solve_result_1, opterror1, solve_result_2, opterror2, solve_result_3, opterror3
end
