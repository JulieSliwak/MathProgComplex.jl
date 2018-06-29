
function SDPInstance_cplx2real(sdp::SDPInstance{T}) where T<:Complex

    block_to_vartype = Dict{String, Symbol}()
    sdpblocks = Dict{Tuple{Moment, String, Exponent, Exponent}, Float64}()
    sdplinsym = Dict{Tuple{Moment, String, Exponent}, Float64}()
    sdplin = Dict{Tuple{Moment, Exponent}, Float64}()
    sdpcst = Dict{Moment, Float64}()

    matrix_terms = Dict{String, Set{Tuple{Exponent, Exponent}}}()

    ## Complex blocks to real
    for ((moment, block_name, γ, δ), coeff) in sdp.blocks
        ctr_re, ctr_im = cplx2real_sdpctr(moment)
        γ_re, γ_im = cplx2real_sdpctr(γ)
        δ_re, δ_im = cplx2real_sdpctr(δ)

        @assert ctr_re != ctr_im

        println("$moment  is $ctr_re, $ctr_im")

        !haskey(matrix_terms, block_name) && (matrix_terms[block_name] = Set{Tuple{Exponent, Exponent}}())
        push!(matrix_terms[block_name], (min(γ, δ), max(γ, δ)))

        # Convert complex linear term to real ones
        ## Real part
        sdpblocks[(ctr_re, block_name, γ_re, δ_re)] = real(coeff)
        sdpblocks[(ctr_re, block_name, γ_im, δ_re)] = -imag(coeff)

        ## Imag part
        if product(moment.conj_part, moment.expl_part) != Exponent()

            sdpblocks[(ctr_im, block_name, γ_re, δ_re)] = imag(coeff)
            sdpblocks[(ctr_im, block_name, γ_im, δ_re)] = real(coeff)
        end
    end

    for (block_name, coords) in matrix_terms
        warn(block_name)
        for (γ, δ) in coords
            @show γ, δ
            γ_re, γ_im = cplx2real_sdpctr(γ)
            δ_re, δ_im = cplx2real_sdpctr(δ)

            Xi_ctrname_re = get_Xictrname_re(block_name, γ, δ)

            sdpblocks[(Xi_ctrname_re, block_name, γ_re, δ_re)] =  1
            sdpblocks[(Xi_ctrname_re, block_name, γ_im, δ_im)] = -1
            info("sdpblocks[($Xi_ctrname_re, $block_name, $γ_re, $δ_re)] =  1")
            info("sdpblocks[($Xi_ctrname_re, $block_name, $γ_im, $δ_im)] = -1")


            Xi_ctrname_im = get_Xictrname_im(block_name, γ, δ)

            sdpblocks[(Xi_ctrname_im, block_name, γ_im, δ_re)] = 1
            sdpblocks[(Xi_ctrname_im, block_name, δ_im, γ_re)] = 1
            info("sdpblocks[($Xi_ctrname_im, $block_name, $γ_im, $δ_re)] = 1")
            info("sdpblocks[($Xi_ctrname_im, $block_name, $δ_im, $γ_re)] = 1")

        end
    end

    ## Complex symetric blocks to real
    for ((moment, block_name, var), coeff) in sdp.linsym
        ctr_re, ctr_im = cplx2real_sdpctr(moment)
        var_re, var_im = cplx2real_sdpctr(var)

        sdplinsym[(ctr_re, block_name, var_re)] = real(coeff)
        sdplinsym[(ctr_re, block_name, var_im)] = -imag(coeff)

        sdplinsym[(ctr_im, block_name, var_re)] = imag(coeff)
        sdplinsym[(ctr_im, block_name, var_im)] = real(coeff)
    end

    ## Complex vectors to real
    for ((moment, var), coeff) in sdp.lin
        ctr_re, ctr_im = cplx2real_sdpctr(moment)
        var_re, var_im = cplx2real_sdpctr(var)

        sdplin[(ctr_re, var_re)] = real(coeff)
        sdplin[(ctr_re, var_im)] = -imag(coeff)

        sdplin[(ctr_im, var_re)] = imag(coeff)
        sdplin[(ctr_im, var_im)] = real(coeff)
    end

    ## Complex coeff to real
    for (moment, coeff) in sdp.cst
        ctr_re, ctr_im = cplx2real_sdpctr(moment)

        sdpcst[ctr_re] = real(coeff)
        sdpcst[ctr_im] = imag(coeff)
    end

    ## Convert variable types
    for (block, vartype) in sdp.block_to_vartype
        if vartype == :SDPC
            block_to_vartype[block] = :SDP
        elseif vartype == :SymC
            block_to_vartype[block] = :Sym
        else
            error("SDPInstance_cplx2real(): Unhandled matrix type $vartype for $block")
        end
    end

    return SDPInstance(block_to_vartype, sdpblocks, sdplinsym, sdplin, sdpcst)
end



function cplx2real_sdpctr(moment::Moment)
    ctr_re = Moment(moment.conj_part, moment.expl_part, moment.clique*"_Re")
    ctr_im = Moment(moment.conj_part, moment.expl_part, moment.clique*"_Im")

    return ctr_re, ctr_im
end

function cplx2real_sdpctr(expo::Exponent)
    expo_expl = expo

    (expo.degree.conjvar > 0) && (expo_expl = conj(expo))
    @assert expo_expl.degree.conjvar == 0

    expo_re, expo_im = Exponent(), Exponent()
    expo_re = Exponent(Variable(string(expo_expl, "_Re"), Real))
    expo_im = Exponent(Variable(string(expo_expl, "_Im"), Real))
    return expo_re, expo_im
end

function get_Xictrname_re(block_name::String, γ::Exponent, δ::Exponent)
    δ_ctr = δ
    if product(γ, δ) == Exponent()
        δ_ctr = Exponent(Variable("1_ctr", Real))
    end
    @show Moment(γ, δ_ctr, block_name*"_ReCtr")
    return Moment(γ, δ_ctr, block_name*"_ReCtr")
    # return Exponent(Variable(block_name*"_Re", Real)), Exponent(Variable(string(γ, "_", δ), Real))
end

function get_Xictrname_im(block_name::String, γ::Exponent, δ::Exponent)
    δ_ctr = δ
    if product(γ, δ) == Exponent()
        δ_ctr = Exponent(Variable("1_ctr", Real))
    end
    @show Moment(γ, δ_ctr, block_name*"_ImCtr")
    return Moment(γ, δ_ctr, block_name*"_ImCtr")
    # return Exponent(Variable(block_name*"_Im", Real)), Exponent(Variable(string(γ, "_", δ), Real))
end