export add!, add, +, -, merge

## Polynomial addition
function add!(p::Polynomial, p1::T) where T<:AbstractPolynomial
    return add!(p, convert(Polynomial, p1))
end

function add!(p::Polynomial, p1::Polynomial)
    for (cur_expo, λ) in p1
        λ != 0 || continue
        add_to_dict!(p.poly, cur_expo, λ)
    end
    p.degree.explvar = max(p.degree.explvar, p1.degree.explvar)
    p.degree.conjvar = max(p.degree.conjvar, p1.degree.conjvar)
end

function add(p1::Polynomial, p2::Polynomial)
    seek_efficiency() && warn("Unefficient implementation\n", @__FILE__, " ", @__LINE__)
    p = deepcopy(p1)
    add!(p, p2)
    return p
end

function +(p1::T, p2::U) where T<:AbstractPolynomial where U<:AbstractPolynomial
    return add(convert(Polynomial, p1), convert(Polynomial, p2))
end
function +(p1::Number, p2::T) where T<:AbstractPolynomial
    return add(convert(Polynomial, p1), convert(Polynomial, p2))
end
function +(p1::T, p2::Number) where T<:AbstractPolynomial
    return add(convert(Polynomial, p1), convert(Polynomial, p2))
end

### Soustraction
function -(p::T) where T<:AbstractPolynomial
    return -1 * convert(Polynomial, p)
end

function -(p1::T, p2::U) where T<:AbstractPolynomial where U<:AbstractPolynomial
    return p1 + (-p2)
end
function -(p1::Number, p2::T) where T<:AbstractPolynomial
    return p1 + (-p2)
end
function -(p1::T, p2::Number) where T<:AbstractPolynomial
    return p1 + (-p2)
end


## Point addition
"""
  add_coord!(pt::Point, var::Variable, val::Number)

  Sparsely add `val` to the `var` coordinate of `pt`.
"""
function add_coord!(pt::Point, var::Variable, val::Number)
  return add_to_dict!(pt.coords, var, val, isdense=pt.isdense)
end


function add!(pt1::Point, pt2::Point)
    for (var, val) in pt2.coords
        add_coord!(pt1, var, val)
    end
end

function add(pt1::Point, pt2::Point)
    seek_efficiency() && warn("Unefficient implementation\n", @__FILE__, " ", @__LINE__)
    pt_coord = deepcopy(pt1.coords)
    pt = Point(pt_coord)
    add!(pt, pt2)
    return pt
end

+(pt1::Point, pt2::Point) = add(pt1, pt2)
-(pt1::Point, pt2::Point) = pt1 + (-1)*pt2

function merge(x::Point ...)
    pt_res = Point()
    for pt in x
        add!(pt_res, pt)
    end
    return pt_res
end