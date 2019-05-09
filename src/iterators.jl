export iterate, length, haskey, keys, values, getindex

## Overloads for transparent iteration on `Exponent`, `Polyomial`, `Point`

## Exponent iterator
iterate(expo::Exponent) = iterate(expo.expo)
iterate(expo::Exponent, state) = iterate(expo.expo, state)
length(expo::Exponent) = length(expo.expo)
haskey(expo::Exponent, key) = haskey(expo.expo, key)
keys(expo::Exponent) = keys(expo.expo)
values(expo::Exponent) = values(expo.expo)
getindex(expo::Exponent, var::Variable) = expo.expo[var]


## Polynomial
iterate(poly::Polynomial) = iterate(poly.poly)
iterate(poly::Polynomial, state) = iterate(poly.poly, state)
length(poly::Polynomial) = length(poly.poly)
haskey(poly::Polynomial, key) = haskey(poly.poly, key)
keys(poly::Polynomial) = keys(poly.poly)
values(poly::Polynomial) = values(poly.poly)
getindex(poly::Polynomial, expo::Exponent) = poly.poly[expo]


## Point
iterate(pt::Point) = iterate(pt.coords)
iterate(pt::Point, state) = iterate(pt.coords, state)
length(pt::Point) = length(pt.coords)
haskey(pt::Point, key) = haskey(pt.coords, key)
keys(pt::Point) = keys(pt.coords)
values(pt::Point) = values(pt.coords)
getindex(pt::Point, var::Variable) = pt.coords[var]
