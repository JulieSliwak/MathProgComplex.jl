export ==, !=, hash, isless

#############################
## Degree
#############################
function isless(deg1::Degree, deg2::Degree)
  return (deg1.explvar<deg2.explvar) || (deg1.explvar==deg2.explvar && deg1.conjvar<deg2.conjvar)
end


#############################
## Variable
#############################
function isless(a::Variable, b::Variable)
  return isless(a.name, b.name)
end


#############################
## Exponent
#############################
function isless(exp1::Exponent, exp2::Exponent)
    exp1_deg = sum(get_sumdegs(exp1))
    exp2_deg = sum(get_sumdegs(exp2))

    # First order level: sum of degrees
    if exp1_deg < exp2_deg
        return true
    elseif exp1_deg > exp2_deg
        return false
    end

    state1 = start(exp1)
    state2 = start(exp2)

    while !done(exp1, state1) && !done(exp2, state2)
        (i1, state1) = next(exp1, state1)
        (i2, state2) = next(exp2, state2)

        if isequal(i1, i2)
            continue
        else
            return isless(i1,i2)
        end
        if done(exp1, state1) && done(exp2, state2)
            return false
        elseif done(exp1, state1) && !done(exp2, state2)
            return true
        elseif !done(exp1, state1) && done(exp2, state2)
            return false
        end
    end

    return false
end


# function isless(exp1::Exponent, exp2::Exponent)
#     # println("\n\n*** isless($exp1, $exp2)")
#     exp1_explsum, exp1_conjsum = get_sumdegs(exp1)
#     exp2_explsum, exp2_conjsum = get_sumdegs(exp2)
#     # exp1_explsum==0 || exp1_conjsum==0 || warn("isless(::Exponent, Exponent): exp1 has expl and conj vars, order may be ill defined...") # TODO
#     # exp2_explsum==0 || exp2_conjsum==0 || warn("isless(::Exponent, Exponent): exp2 has expl and conj vars, order may be ill defined...")
#     exp1_deg = exp1_explsum + exp1_conjsum
#     exp2_deg = exp2_explsum + exp2_conjsum

#     if exp1_deg < exp2_deg
#         # println("****return true")
#         return true
#     elseif exp1_deg == exp2_deg
#         vars = SortedSet(keys(exp1))
#         union!(vars, keys(exp2))
#         # @show vars
#         for var in vars
#             # @show var
#             # println("->here!")
#             if !haskey(exp2, var) # hence haskey(exp1, var)
#                 # println("->here!1")
#                 # println("****return true")
#                 return true
#             elseif !haskey(exp1, var) # hence haskey(exp2, var)
#                 # println("->here!2")
#                 # println("****return false")
#                 return false
#             elseif exp1[var] < exp2[var]
#                 # println("->here!3")
#                 # println("****return true")
#                 return true
#             elseif exp1[var] == exp2[var]
#                 # println("->here!4")
#                 continue
#             else
#                 # println("->here!5")
#                 # println("****return false")
#                 return false
#             end
#         end
#     else
#         # println("****return false")
#         return false
#     end
#     # println("****return false")
#     return false
# end

#############################
## Polynomial
#############################
# NOTE: order on polynomials ?
