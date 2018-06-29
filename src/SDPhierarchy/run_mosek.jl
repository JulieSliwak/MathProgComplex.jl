using Mosek

printstream(msg::String) = print(msg)

function get_SDPtriplets(problem::SDP_Problem; debug = false)

  nzc = 0
  nza = 0
  for ((objctr, block, var1, var2), coeff) in problem.matrices

      if objctr in problem.obj_keys
        nzc += 1
      else
        nza += 1
      end
  end
  # println("nza : ", nza)
  # println("nzc : ", nzc)
  barai = Int32[0 for i in 1:nza]
  baraj = Int32[0 for i in 1:nza]
  barak = Int32[0 for i in 1:nza]
  baral = Int32[0 for i in 1:nza]
  baraijkl = Float64[0 for i in 1:nza]

  barcj = Int32[0 for i in 1:nzc]
  barck = Int32[0 for i in 1:nzc]
  barcl = Int32[0 for i in 1:nzc]
  barcjkl = Float64[0 for i in 1:nzc]

  nzc=0
  nza=0
  for ((objctr, block, var1, var2), coeff) in problem.matrices
    sdp_block = problem.name_to_sdpblock[block]
    lower = min(sdp_block.var_to_id[var1], sdp_block.var_to_id[var2])
    upper = max(sdp_block.var_to_id[var1], sdp_block.var_to_id[var2])
    if objctr in problem.obj_keys
      nzc+=1
      barcj[nzc] = sdp_block.id
      barck[nzc] = upper
      barcl[nzc] = lower
      barcjkl[nzc] = coeff
    else
      nza+=1
      barai[nza] = problem.name_to_ctr[objctr][1]
      baraj[nza] = sdp_block.id
      barak[nza] = upper
      baral[nza] = lower
      baraijkl[nza] = coeff
    end
  end

  if debug
    println("*********************************************************************************")
    println("Debug -> Reading SDP_Problem")
    @printf("%5s  %5s  %5s  %s\n", "barcj", "barck", "barcl", "barcjkl")
    for i=1:length(barcj)
        @printf("%5i  %5i  %5i  %f\n", barcj[i], barck[i], barcl[i], barcjkl[i])
    end

    @printf("%5s  %5s  %5s  %5s  %s\n", "barai", "baraj", "barak", "baral", "baraijkl")
    for i=1:length(barai)
        @printf("%5i  %5i  %5i  %5i  %f\n", barai[i], baraj[i], barak[i], baral[i], baraijkl[i])
    end
  end
  return barcj, barck, barcl, barcjkl, barai, baraj, barak, baral, baraijkl
end

function get_linterms(problem; debug=debug)
  nza, nzc = 0, 0
  for (objctr, var) in keys(problem.linear)
    if objctr in problem.obj_keys
      nzc += 1
    else
      nza += 1
    end
  end

  ai = zeros(Int32, nza)
  aj = zeros(Int32, nza)
  aij = zeros(Float64, nza)
  cj = zeros(Int32, nzc)
  cjval = zeros(Float64, nzc)

  nza, nzc = 0, 0
  for ((objctr, var), coeff) in problem.linear
    if objctr in problem.obj_keys
      nzc += 1
      cj[nzc] = problem.scalvar_to_id[var]
      cjval[nzc] = coeff
    else
      nza += 1
      ai[nza] = problem.name_to_ctr[objctr][1]
      aj[nza] = problem.scalvar_to_id[var]
      aij[nza] = coeff
    end
  end

  if debug
    warn("--- cj :")
    for i=1:length(cj)
      @printf("%i  %f\n", cj[i], cjval[i])
    end
    warn("--- a :")
    for i=1:length(aj)
      @printf("%i  %i  %f\n", ai[i], aj[i], aij[i])
    end
  end
  return cj, cjval, ai, aj, aij
end

function get_ctrbounds(problem::SDP_Problem; debug = false)
  numcon=length(problem.name_to_ctr)
  MOSEK_KIND = Dict(["EQ"=>MSK_BK_FX, "GEQ"=>MSK_BK_LO, "LEQ"=>MSK_BK_UP, "RNG"=>MSK_BK_RA])
  bkc = Boundkey[ Mosek.Boundkey(1)  for kv in problem.name_to_ctr]
  # constraint
  buc = Float64[0 for i in 1:numcon]
  blc = Float64[0 for i in 1:numcon]
  for (ctrname, ctr) in problem.name_to_ctr
    # @printf("%10d%20s%10d%20s\n", ctr[1], ctrname, MOSEK_KIND[ctr[2]], ctr[2])
    id_ctr=ctr[1]
    lb = ctr[3]
    ub = ctr[4]
    cst = 0
    if haskey(problem.cst_ctr, ctrname)
      cst = problem.cst_ctr[ctrname]
    end
    bkc[id_ctr] = MOSEK_KIND[ctr[2]]
    if bkc[id_ctr] == MSK_BK_UP
      buc[id_ctr] = ub - cst
    elseif bkc[id_ctr] == MSK_BK_LO
      blc[id_ctr] = lb - cst
    elseif bkc[id_ctr] == MSK_BK_FX
      blc[id_ctr] = lb - cst
      buc[id_ctr] = ub - cst
    elseif bkc[id_ctr] == MSK_BK_RA
      blc[id_ctr] = lb - cst
      buc[id_ctr] = ub - cst
    else
      error("get_ctrbounds() : Unknown constraint kind $(ctr[2]) $(bkc[id_ctr]) $(MSK_BK_FX[1])")
    end
  end
  if debug
    warn("get_ctrbounds(): done")
    @show numcon
    @show bkc
    @show blc
    @show buc
  end
  return numcon, bkc, blc, buc
end

function get_varbounds(problem::SDP_Problem)
  MSK_INFINITY = 1.0e30

  numvar = length(problem.scalvar_to_id)

  sub = [i for i in 1:numvar]
  bkx = [MSK_BK_FR for i in 1:numvar]
  blx = [-MSK_INFINITY for i in 1:numvar]
  bux = [MSK_INFINITY for i in 1:numvar]

  return sub, bkx, blx, bux
end

function solve_mosek(problem::SDP_Problem, primal::SortedDict{Tuple{String,String,String}, Float64},
                                           dual::SortedDict{Tuple{String, String, String}, Float64};
                                           debug = false,
                                           logname = "",
                                           printlog = true,
                                           msk_maxtime = -1,            #Default -1 means no time limit
                                           sol_info = OrderedDict())
  empty!(primal)
  empty!(dual)
  primobj = NaN
  dualobj = NaN

  # Get nb of variables, dimensions, nb of ctrs, bounds.
  nbarvar = length(problem.id_to_sdpblock)
  barvardim = [ length(problem.id_to_sdpblock[block].var_to_id) for block in 1:nbarvar ]

  numvar = length(problem.scalvar_to_id)

  numcon, bkc, blc, buc = get_ctrbounds(problem)
  sub, bkx, blx, bux = get_varbounds(problem)

  # Get objective and ctr matrix contribution
  barcj, barck, barcl, barcjkl, barai, baraj, barak, baral, baraijkl = get_SDPtriplets(problem, debug=debug)

  # Get obj and ctr linear contribution
  cj, cjval, ai, aj, aij = get_linterms(problem, debug=debug)

  # Create a task object and attach log stream printer
  maketask() do task
      logname != "" && linkfiletostream(task, MSK_STREAM_LOG, logname, 0)

      printlog && putstreamfunc(task,MSK_STREAM_LOG,printstream)

      # Append SDP matrix variables and scalar variables.
      # The variables will initially be fixed at zero.
      appendbarvars(task,barvardim)
      appendvars(task, numvar)

      putvarboundlist(task, sub, bkx, blx, bux)

      # Append 'numcon' empty constraints.
      # The constraints will initially have no bounds.
      appendcons(task,numcon)
      # Set the bounds on constraints.
      putconboundslice(task,1,numcon+1, bkc,blc,buc)

      # Minimize
      putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

      # Set constraints SDP vars coeffs
      putbarablocktriplet(task, length(barai), barai, baraj, barak, baral, baraijkl)

      # Set contraints Sym vars coeffs as linear part
      putaijlist(task, ai, aj, aij)

      # Objective matrices and constant
      putbarcblocktriplet(task, length(barcj), barcj, barck, barcl, barcjkl)
      putclist(task, cj, cjval)

      obj_cst = 0
      for objmoment in problem.obj_keys
        if haskey(problem.cst_ctr, objmoment)
          obj_cst += problem.cst_ctr[objmoment]
        end
      end
      obj_cst != 0 && putcfix(task, obj_cst)

      # putintparam(task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_NONE)
      # putintparam(task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_AGGRESSIVE)

      # putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_DFEAS, 1e-12)
      # putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_INFEAS, 1e-12)
      # putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_MU_RED, 1e-12)
      # putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_PFEAS, 1e-12)
      # putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, 1e-1)

      # Set Mosek maxtime, default is infinity
      putdouparam(task, MSK_DPAR_OPTIMIZER_MAX_TIME, msk_maxtime)

      if debug
        dump_mosek_model(task)
      end

      optimize(task)
      solutionsummary(task,MSK_STREAM_MSG)

      # Get status information about the solution
      prosta = getprosta(task, MSK_SOL_ITR)
      solsta = getsolsta(task, MSK_SOL_ITR)

      save_SDPsolinfo!(task, sol_info)

      if solsta == MSK_SOL_STA_OPTIMAL || solsta == MSK_SOL_STA_NEAR_OPTIMAL
        # Get objective
        primobj = getprimalobj(task, MSK_SOL_ITR)
        dualobj = getdualobj(task, MSK_SOL_ITR)

        # Get primal and dual solutions
        get_primalsol!(task, problem, primal)
        get_dualsol!(task, problem, dual)


      elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
          println("Dual infeasibility.\n")
      elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
          println("Primal infeasibility.\n")
      elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
          println("Near dual infeasibility.\n")
      elseif  solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
          println("Near primal infeasibility.\n")
      elseif  solsta == MSK_SOL_STA_UNKNOWN
          println("Unknown solution status.\n")
      else
          println("Other solution status.\n")
      end
  end

  return primobj, dualobj
end


function save_SDPsolinfo!(task::Mosek.Task, sol_info::OrderedDict)
  prosta = getprosta(task, MSK_SOL_ITR)
  solsta = getsolsta(task, MSK_SOL_ITR)

  sol_info[:slv_prosta] = get_prosta_map()[prosta]
  sol_info[:slv_solsta] = get_solsta_map()[solsta]

  if solsta == MSK_SOL_STA_OPTIMAL || solsta == MSK_SOL_STA_NEAR_OPTIMAL
    sol_info[:slv_primobj] = getprimalobj(task, MSK_SOL_ITR)
    sol_info[:slv_dualobj] = getdualobj(task, MSK_SOL_ITR)
    sol_info[:slv_primfeas] = getdouinf(task, MSK_DINF_INTPNT_PRIMAL_FEAS)
    sol_info[:slv_dualfeas] = getdouinf(task, MSK_DINF_INTPNT_DUAL_FEAS)
    sol_info[:slv_solvetime] = getdouinf(task, MSK_DINF_OPTIMIZER_TIME)
  else
    sol_info[:slv_primobj] = -1.0
    sol_info[:slv_dualobj] = -1.0
    sol_info[:slv_primfeas] = -1.0
    sol_info[:slv_dualfeas] = -1.0
    sol_info[:slv_solvetime] = -1.0
  end
end


function get_primalsol!(task::Mosek.Task, problem::SDP_Problem,
                                          primal::SortedDict{Tuple{String,String,String}, Float64})
  for (id, block) in problem.id_to_sdpblock
    barx = getbarxj(task, MSK_SOL_ITR, id)
    all_variables = ["" for kv in block.var_to_id]
    for (var, varid) in block.var_to_id
      all_variables[varid] = var
    end
    n = 0
    for j in 1:length(all_variables)
      for i in j:length(all_variables)
        n+=1
        # @printf("%15s%15s%15s%25.10f\n", id_block[2].name, all_variables[i], all_variables[j], barx[n])
        primal[block.name, all_variables[i], all_variables[j]] = barx[n]
      end
    end
  end
end


function get_dualsol!(task::Mosek.Task, problem::SDP_Problem,
                                        dual::SortedDict{Tuple{String, String, String}, Float64})
  for (blockid, block) in problem.id_to_sdpblock
    bars = getbarsj(task, MSK_SOL_ITR, blockid)
    id_to_var = OrderedDict([id=>var for (var,id) in block.var_to_id])
    it = 0
    for j=1:length(id_to_var), i=j:length(id_to_var)
      it += 1
      # println("($(block.name), $(id_to_var[i]), $(id_to_var[j])) = $(bars[it])")
      dual[(block.name, id_to_var[i], id_to_var[j])] = bars[it]
    end
  end
end


function dump_mosek_model(task)
  println("*********************************************************************************")
  println("Debug -> Reading Mosek params")
  println("MSK_DPAR_INTPNT_CO_TOL_DFEAS,  $(getdouparam(task, MSK_DPAR_INTPNT_CO_TOL_DFEAS))")
  println("MSK_DPAR_INTPNT_CO_TOL_INFEAS,   $(getdouparam(task, MSK_DPAR_INTPNT_CO_TOL_INFEAS))")
  println("MSK_DPAR_INTPNT_CO_TOL_MU_RED,   $(getdouparam(task, MSK_DPAR_INTPNT_CO_TOL_MU_RED))")
  println("MSK_DPAR_INTPNT_CO_TOL_PFEAS,  $(getdouparam(task, MSK_DPAR_INTPNT_CO_TOL_PFEAS))")
  println("MSK_DPAR_INTPNT_CO_TOL_REL_GAP,  $(getdouparam(task, MSK_DPAR_INTPNT_CO_TOL_REL_GAP))")

  println("\nDebug -> Reading Mosek problem")
  num, subcj, subck, subcl, valcjkl = getbarcblocktriplet(task)
  @printf("%5s  %5s  %5s  %s\n", "subcj", "subck", "subcl", "valcjkl")
  for i=1:num
      @printf("%5i  %5i  %5i  %f\n", subcj[i], subck[i], subcl[i], valcjkl[i])
  end

  num, subai, subaj, subak, subal, valajkl = getbarablocktriplet(task)
  @printf("\n%5s  %5s  %5s  %5s  %s\n", "subai", "subaj", "subak", "subal", "valajkl")
  for i=1:num
      @printf("%5i  %5i  %5i  %5i  %f\n", subai[i], subaj[i], subak[i], subal[i], valajkl[i])
  end

  boundkeys, lbs, ubs = getconboundslice(task, 1, numcon+1)
  @show boundkeys
  @show lbs
  @show ubs
  println("*********************************************************************************")
end

function get_solsta_map()
  return OrderedDict(MSK_SOL_STA_UNKNOWN => "Status of the solution is unknown.",
                     MSK_SOL_STA_OPTIMAL => "The solution is optimal.",
                     MSK_SOL_STA_PRIM_FEAS => "The solution is primal feasible.",
                     MSK_SOL_STA_DUAL_FEAS => "The solution is dual feasible.",
                     MSK_SOL_STA_PRIM_AND_DUAL_FEAS => "The solution is both primal and dual feasible.",
                     MSK_SOL_STA_NEAR_OPTIMAL => "The solution is nearly optimal.",
                     MSK_SOL_STA_NEAR_PRIM_FEAS => "The solution is nearly primal feasible.",
                     MSK_SOL_STA_NEAR_DUAL_FEAS => "The solution is nearly dual feasible.",
                     MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS => "The solution is nearly both primal and dual feasible.",
                     MSK_SOL_STA_PRIM_INFEAS_CER => "The solution is a certificate of primal infeasibility.",
                     MSK_SOL_STA_DUAL_INFEAS_CER => "The solution is a certificate of dual infeasibility.",
                     MSK_SOL_STA_NEAR_PRIM_INFEAS_CER => "The solution is almost a certificate of primal infeasibility.",
                     MSK_SOL_STA_NEAR_DUAL_INFEAS_CER => "The solution is almost a certificate of dual infeasibility.",
                     MSK_SOL_STA_PRIM_ILLPOSED_CER => "The solution is a certificate that the primal problem is illposed.",
                     MSK_SOL_STA_DUAL_ILLPOSED_CER => "The solution is a certificate that the dual problem is illposed.",
                     MSK_SOL_STA_INTEGER_OPTIMAL => "The primal solution is integer optimal.",
                     MSK_SOL_STA_NEAR_INTEGER_OPTIMAL => "The primal solution is near integer optimal.")
end

function get_prosta_map()
  return OrderedDict(MSK_PRO_STA_UNKNOWN => "Unknown problem status.",
                     MSK_PRO_STA_PRIM_AND_DUAL_FEAS => "The problem is primal and dual feasible.",
                     MSK_PRO_STA_PRIM_FEAS => "The problem is primal feasible.",
                     MSK_PRO_STA_DUAL_FEAS => "The problem is dual feasible.",
                     MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS => "The problem is at least nearly primal and dual feasible.",
                     MSK_PRO_STA_NEAR_PRIM_FEAS => "The problem is at least nearly primal feasible.",
                     MSK_PRO_STA_NEAR_DUAL_FEAS => "The problem is at least nearly dual feasible.",
                     MSK_PRO_STA_PRIM_INFEAS => "The problem is primal infeasible.",
                     MSK_PRO_STA_DUAL_INFEAS => "The problem is dual infeasible.",
                     MSK_PRO_STA_PRIM_AND_DUAL_INFEAS => "The problem is primal and dual infeasible.",
                     MSK_PRO_STA_ILL_POSED => "The problem is ill-posed. For example, it may be primal and dual feasible but have a positive duality gap.",
                     MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED => "The problem is either primal infeasible or unbounded. This may occur for mixed-integer problems.")
end
