# export export_SDP_Instance

function export_SDP_Instance(sdp_pb::SDP_Problem, path::String)
    # Export blocks of constraints
    blocks_file = joinpath(path, "matrix.sdp")
    !isfile(blocks_file) || rm(blocks_file)

    open(blocks_file, "a") do fblocks
        print_matrices(fblocks, sdp_pb.matrices)
    end

    # Export linear part of constraints
    lin_file = joinpath(path, "lin.sdp")
    !isfile(lin_file) || rm(lin_file)

    open(lin_file, "a") do flin
        print_scalarlinear(flin, sdp_pb.linear)
    end

    # Export constants of constraints
    cst_file = joinpath(path, "const.sdp")
    !isfile(cst_file) || rm(cst_file)

    open(cst_file, "a") do fcst
        print_constant(fcst, sdp_pb.cst_ctr)
    end

    # Export bloc types
    types_file = joinpath(path, "types.sdp")
    !isfile(types_file) || rm(types_file)

    open(types_file, "a") do ftypes
        print_sdpblocks(ftypes, sdp_pb.name_to_sdpblock)
    end
end

function print_matrices(io::IO, matrices::SortedDict{Tuple{SDP_Moment, String, String, String}, Float64}, indentedprint=true)
    MPC.print_blocksfile_header(io)

    ctrkey1len =  max( maximum(x->length(x[1][1]), keys(matrices)), length("#j_conj"))
    ctrkey2len =  max( maximum(x->length(x[1][2]), keys(matrices)), length("j_expl"))
    ctrkey3len = max( maximum(x->length(x[1][3]), keys(matrices)), length("clique"))
    blocklen =  max( maximum(x->length(x[2]), keys(matrices)), length("Zi"))
    rowlen =    max( maximum(x->length(x[3]), keys(matrices)), length("k"))
    collen =    max( maximum(x->length(x[4]), keys(matrices)), length("l"))
    MPC.print_string(io, "#j_conj", ctrkey1len, indentedprint=indentedprint)
    MPC.print_string(io, "j_expl", ctrkey2len, indentedprint=indentedprint)
    MPC.print_string(io, "clique", ctrkey3len, indentedprint=indentedprint)
    MPC.print_string(io, "Zi", blocklen, indentedprint=indentedprint)
    MPC.print_string(io, "k", rowlen, indentedprint=indentedprint)
    MPC.print_string(io, "l", collen, indentedprint=indentedprint)
    MPC.print_string(io, "Real(A_ij[k, l])", 23, indentedprint=indentedprint)
    MPC.print_string(io, "Imag(A_ij[k, l])", 23, indentedprint=indentedprint)
    println(io)

    for (((ctrname1, ctrname2, ctrname3), blockname, var1, var2), λ) in matrices

        MPC.print_string(io, ctrname1, ctrkey1len, indentedprint=indentedprint)
        MPC.print_string(io, ctrname2, ctrkey2len, indentedprint=indentedprint)
        MPC.print_string(io, ctrname3, ctrkey3len, indentedprint=indentedprint)
        MPC.print_string(io, blockname, blocklen, indentedprint=indentedprint)
        MPC.print_string(io, var1, rowlen, indentedprint=indentedprint)
        MPC.print_string(io, var2, collen, indentedprint=indentedprint)
        @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
    end
end

function print_scalarlinear(io::IO, linear::SortedDict{Tuple{SDP_Moment, String}, Float64}, indentedprint=true)
    MPC.print_linfile_header(io)

    (length(linear) == 0) && return

    ctrkey1len =  max( maximum(x->length(x[1][1]), keys(linear)), length("#j_conj"))
    ctrkey2len =  max( maximum(x->length(x[1][2]), keys(linear)), length("j_expl"))
    ctrkey3len = max( maximum(x->length(x[1][3]), keys(linear)), length("clique"))
    varlen =    max( maximum(x->length(x[2]), keys(linear)), length("x[k]"))
    MPC.print_string(io, "#j_conj", ctrkey1len, indentedprint=indentedprint)
    MPC.print_string(io, "j_expl", ctrkey2len, indentedprint=indentedprint)
    MPC.print_string(io, "clique", ctrkey3len, indentedprint=indentedprint)
    MPC.print_string(io, "k", varlen, indentedprint=indentedprint)
    MPC.print_string(io, "Real(b_j[k])", 23, indentedprint=indentedprint)
    MPC.print_string(io, "Imag(b_j[k])", 23, indentedprint=indentedprint)
    println(io)

    for (((ctrname1, ctrname2, ctrname3), varname), λ) in linear

        MPC.print_string(io, ctrname1, ctrkey1len, indentedprint=indentedprint)
        MPC.print_string(io, ctrname2, ctrkey2len, indentedprint=indentedprint)
        MPC.print_string(io, ctrname3, ctrkey3len, indentedprint=indentedprint)
        MPC.print_string(io, varname, varlen, indentedprint=indentedprint)
        @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
    end
end

function print_constant(io::IO, cst::SortedDict{SDP_Moment, Float64}, indentedprint=true)
    MPC.print_cstfile_header(io)

    (length(cst) == 0) && return

    ctrkey1len =  max( maximum(x->length(x[1]), keys(cst)), length("#j_conj"))
    ctrkey2len =  max( maximum(x->length(x[2]), keys(cst)), length("j_expl"))
    ctrkey3len = max( maximum(x->length(x[3]), keys(cst)), length("clique"))

    MPC.print_string(io, "#j_conj", ctrkey1len, indentedprint=indentedprint)
    MPC.print_string(io, "j_expl", ctrkey2len, indentedprint=indentedprint)
    MPC.print_string(io, "clique", ctrkey3len, indentedprint=indentedprint)
    MPC.print_string(io, "Real(c_j)", 23, indentedprint=indentedprint)
    MPC.print_string(io, "Imag(c_j)", 23, indentedprint=indentedprint)
    println(io)

    for ((ctrname1, ctrname2, ctrname3), λ) in cst

        MPC.print_string(io, ctrname1, ctrkey1len, indentedprint=indentedprint)
        MPC.print_string(io, ctrname2, ctrkey2len, indentedprint=indentedprint)
        MPC.print_string(io, ctrname3, ctrkey3len, indentedprint=indentedprint)
        @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
    end
end


function print_sdpblocks(io::IO, name_to_sdpblock::SortedDict{String, SDP_Block}, indentedprint=true)
    MPC.print_typesfile_header(io)

    cstrlen = max( maximum(x->length(x), keys(name_to_sdpblock)), length("#Zi"))

    MPC.print_string(io, "#Zi", cstrlen, alignright=false); println(io, " type")

    for (blockname, block) in name_to_sdpblock
        MPC.print_string(io, blockname, cstrlen)
        println(io, " SDP")
    end
end