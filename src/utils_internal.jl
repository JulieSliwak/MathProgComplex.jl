
"""
    add_to_dict!(dict::SortedDict{Any, V}, key, val::V; isdense=false) where V<:Number

    *Sparsely* add `val` to the `key` entry of `dict` dictionnary (if not `isdense`). 
    That is creates the entry if needed, deletes it if the resulting value is null.
"""
function add_to_dict!(dict::SortedDict{U, V}, key::U, val::T; isdense = false) where T<:Number where U where V<:Number
    if !haskey(dict, key)
        dict[key] = convert(V, 0)
    end
    dict[key] += convert(V, val)

    if (dict[key]==0) && !isdense
        delete!(dict, key)
    end
end