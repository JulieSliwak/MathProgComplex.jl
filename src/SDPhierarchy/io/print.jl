export print

###############################################################################
####  Relaxation context
###############################################################################
function print(io::IO, relctx::RelaxationContext)
    print(io, "RelaxationContext:\n")
    print(io, "ismultiordered         : $(relctx.ismultiordered)\n")
    print(io, "issparse               : $(relctx.issparse)\n")
    print(io, "symmetries             : $(relctx.symmetries)\n")
    print(io, "hierarchykind          : $(relctx.hierarchykind)\n")
    print(io, "renamevars             : $(relctx.renamevars)\n")
    for cstrname in sort(collect(keys(relctx.di)))
        di = relctx.di[cstrname]
        print(io, "di                     : $cstrname  \t=> $di\n")
    end
    for cstrname in sort(collect(keys(relctx.ki)))
        ki = relctx.ki[cstrname]
        print(io, "ki                     : $cstrname  \t=> $ki\n")
    end
    for cstrname in sort(collect(keys(relctx.cstrtypes)))
        cstrtype = relctx.cstrtypes[cstrname]
        print(io, "bar var types          : $cstrname  \t=> $(string(cstrtype))\n")
    end
end


###############################################################################
####  SDPDual
###############################################################################
function print(io::IO, momentrelax::SDPDual{T}) where T
    println(io, "Moment Relaxation Problem:")
    println(io, "→ Objective: ")
    momentlen = maximum(x->length(string(x)), keys(momentrelax.objective))
    for moment in sort(collect(keys(momentrelax.objective)))
        coeff = momentrelax.objective[moment]
        print_string(io, string(moment), momentlen)
        println(io, " $coeff")
    end

    println(io, "→ Constraints:")
    for (cstrname, blocname) in sort(collect(keys(momentrelax.constraints)))
        mmtmat = momentrelax.constraints[(cstrname, blocname)]
        println(io, " → $cstrname, $blocname")
        println(io, mmtmat)
    end

    println(io, "→ Moments clique overlap:")
    if length(momentrelax.moments_overlap) > 0
        mmtlength = maximum(x->length(string(x)), keys(momentrelax.moments_overlap))
        for moment in sort(collect(keys(momentrelax.moments_overlap)))
            cliquenames = momentrelax.moments_overlap[moment]
            print(io, " → ")
            print_string(io, string(moment), mmtlength)
            for clique in sort(collect(cliquenames)) print(io, "$clique, ") end
            @printf(io, "\b\b \n")
        end
    else
        print(io, "  None")
    end
end


###############################################################################
####  SDPPrimal
###############################################################################
function print(io::IO, sdpinst::SDPPrimal)
    println(io, " -- SDP Blocks:")
    print(io, sdpinst.blocks)
    println(io, " -- linear part:")
    print(io, sdpinst.lin, sdpinst.linsym)
    println(io, " -- const part:")
    print(io, sdpinst.cst)
    println(io, " -- mat var types:")
    for blockname in sort(collect(keys(sdpinst.block_to_vartype)))
        blocktype = sdpinst.block_to_vartype[blockname]
        println(io, "   $blockname  \t $blocktype")
    end
end


function print(io::IO, sdpblocks::DictType{Tuple{Moment, String, Exponent, Exponent}, T}; indentedprint=true) where T
    print_blocksfile(io, sdpblocks; indentedprint=indentedprint, print_header=false)
end

function print(io::IO, sdplin::DictType{Tuple{Moment, Exponent}, T}, sdplinsym::DictType{Tuple{Moment, String, Exponent}, T}; indentedprint=true) where T
    print_linfile(io, sdplin, sdplinsym; indentedprint=indentedprint, print_header=false)
end


function print(io::IO, sdpcst::DictType{Moment, T}; indentedprint=true) where T
    print_cstfile(io, sdpcst; indentedprint=indentedprint, print_header=false)
end


###############################################################################
####  SDP_Problem
###############################################################################
function print(io::IO, sdp::SDP_Problem)
  for cstr in keys(sdp.name_to_sdpblock)
    block = sdp.name_to_sdpblock[cstr]
    println(io, "  sdp   : $cstr -> $block")
  end

  println(io, "  objk  : $(sdp.obj_keys)")

  for name in keys(sdp.name_to_ctr)
    ctr = sdp.name_to_ctr[name]
    println(io, "  ctr   : $name \t $ctr")
  end

  for name in keys(sdp.matrices)
    mat = sdp.matrices[name]
    println(io, "  matrix: $name \t $mat")
  end

  for name in keys(sdp.linear)
    lin = sdp.linear[name]
    println(io, "  lin   : $name \t $lin")
  end

  for name in keys(sdp.cst_ctr)
    cst = sdp.cst_ctr[name]
    println(io, "  cst   : $name \t $cst")
  end
end