export print

function print(io::IO, sdp::SDP_Problem)
  for cstr in sort(keys(sdp.name_to_sdpblock))
    block = sdp.name_to_sdpblock[cstr]
    println(io, "  sdp   : $cstr -> $block")
  end

  println(io, "  objk  : $(sdp.obj_keys)")

  for name in sort(keys(sdp.name_to_ctr))
    ctr = sdp.name_to_ctr[name]
    println(io, "  ctr   : $name \t $ctr")
  end

  for name in sort(keys(sdp.matrices))
    mat = sdp.matrices[name]
    println(io, "  matrix: $name \t $mat")
  end

  for name in sort(keys(sdp.linear))
    lin = sdp.linear[name]
    println(io, "  lin   : $name \t $lin")
  end

  for name in sort(keys(sdp.cst_ctr))
    cst = sdp.cst_ctr[name]
    println(io, "  cst   : $name \t $cst")
  end
end