
function Moment(expo::Exponent, clique::String)
    α, β = Exponent(), Exponent()

    for (var, deg) in expo
        product!(α, Exponent(SortedDict(var=>Degree(0, deg.conjvar))))
        product!(β, Exponent(SortedDict(var=>Degree(deg.explvar, 0))))
    end
    return Moment(α, β, clique)
end

hash(mom::Moment, h::UInt) = hash(mom.conj_part, hash(mom.expl_part, hash(mom.clique, h)))

isless(mom1::Moment, mom2::Moment) = isless((mom1.conj_part, mom1.expl_part, mom1.clique),
                                            (mom2.conj_part, mom2.expl_part, mom2.clique))

function print(io::IO, mom::Moment)
    print(io, "$(mom.expl_part * mom.conj_part) ($(mom.clique))")
end

==(mom1::Moment, mom2::Moment) = ((mom1.conj_part, mom1.expl_part, mom1.clique) == (mom2.conj_part, mom2.expl_part, mom2.clique))
