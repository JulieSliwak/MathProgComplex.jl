export build_SDP_Instance_from_sdpfiles

function build_SDP_Instance_from_sdpfiles(path::String)
  sdp_read = read_SDPPrimal(path)

  # println("VAR_TYPES size:     $(size(sdp_read.VAR_TYPES))")
  # println("BLOCKS size:        $(size(sdp_read.BLOCKS))")
  # println("LINEAR size:        $(size(sdp_read.LINEAR))")
  # println("CONST size:         $(size(sdp_read.CONST))")

  sdp_inst = SDP_Problem()

  # Build numerical descirption of problem, indexed on strings
  set_matrices!(sdp_inst, sdp_read)
  set_linear!(sdp_inst, sdp_read)
  set_const!(sdp_inst, sdp_read)

  # Build structural information and string to id maps
  set_vartypes!(sdp_inst, sdp_read)

  set_blocks!(sdp_inst)
  set_scalvars!(sdp_inst)


  set_constraints!(sdp_inst)

  return sdp_inst
end



function read_SDPPrimal(path::String)
  BLOCKS = readdlm(joinpath(path, "matrix.sdp"), String)
  if isfile(joinpath(path, "lin.sdp")) && length(matchall(r"\n", readstring(joinpath(path, "lin.sdp")))) > 7
    LINEAR = readdlm(joinpath(path, "lin.sdp"), String)
  else
    LINEAR = []
  end
  CONST = readdlm(joinpath(path, "const.sdp"), String)
  VAR_TYPES = readdlm(joinpath(path, "types.sdp"), String)

  SDP_Instance(VAR_TYPES,
               BLOCKS,
               LINEAR,
               CONST)
end

function set_matrices!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  for i=1:size(instance.BLOCKS, 1)
    ctr_name = (instance.BLOCKS[i, 1], instance.BLOCKS[i, 2], instance.BLOCKS[i, 3])
    (block_name, var1, var2, coeff) = instance.BLOCKS[i, 4:7]

    # Sort variables for triangular matrix storage
    var1, var2 = min(var1, var2), max(var1, var2)

    # if haskey(sdp.name_to_sdpblock, block_name)
    if !haskey(sdp.matrices, (ctr_name, block_name, var1, var2))
      sdp.matrices[(ctr_name, block_name, var1, var2)] = parse(coeff)
    else
      error("set_matrices!(): sdp.matrices already has key ($ctr_name, $block_name, $var1, $var2) with val $(sdp.matrices[(ctr_name, block_name, var1, var2)]), $(parse(coeff))")
    end

    # else
    #   error("set_matrices!(): Unhandled matrix var $block_name")
    # end
  end

  if debug
  end
end


function set_linear!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  if length(instance.LINEAR) != 0
    for i=1:size(instance.LINEAR, 1)
      ctr_name = (instance.LINEAR[i, 1], instance.LINEAR[i, 2], instance.LINEAR[i, 3])
      (var, coeff) = instance.LINEAR[i, 4:5]

      if !haskey(sdp.linear, (ctr_name, var))
        sdp.linear[(ctr_name, var)] = parse(coeff)
      else
        error("set_linear!(): sdp.linear already has key ($ctr_name, $var) with val $(sdp.linear[(ctr_name, var)]). New val is $(parse(coeff))")
      end
    end
  end

  if debug
  end
end

function set_const!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  for i=1:size(instance.CONST, 1)
    ctr_name = (instance.CONST[i, 1], instance.CONST[i, 2], instance.CONST[i, 3])
    coeff = instance.CONST[i, 4]

    @assert !haskey(sdp.cst_ctr, ctr_name)
    if coeff != 0
      sdp.cst_ctr[ctr_name] = parse(coeff)
    end
  end

  if debug
  end
end


"""
set_vartypes!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)

Input all matrix varibales structural information regarding name and kind in the appropriate attributes of `SDP_Problem`.
"""
function set_vartypes!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  n_sdp = 0
  for i=1:size(instance.VAR_TYPES, 1)
    (block_name, block_type) = instance.VAR_TYPES[i, :]

    if block_type == "SDP"
      n_sdp += 1
      block = SDP_Block(n_sdp, block_name)
      sdp.name_to_sdpblock[block_name] = block
      sdp.id_to_sdpblock[n_sdp] = block

    else
      error("set_vartypes()!: Unknown blockvar type $block_type for variable $block_name")
    end
  end

  if debug
  end
end
