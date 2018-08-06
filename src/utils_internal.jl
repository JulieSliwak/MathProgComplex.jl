export seek_efficiency, seek_efficiency!

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

    if !isdense && (dict[key]==0)
        delete!(dict, key)
    end
end

"""
    seek_efficiency()

Get the `seek_efficiency` flag.
If set to true, warnings will be display by inefficient functions: use of
`deepcopy`, `add` instead of `add!`...
"""
seek_efficiency() = SEEK_EFFICIENCY[1]

"""
    seek_efficiency!(arg)

Set the `seek_efficiency` flag.
If set to true, warnings will be display by inefficient functions: use of
`deepcopy`, `add` instead of `add!`...
"""
function seek_efficiency!(arg::Bool)
    SEEK_EFFICIENCY[1] = arg
end