export start, next, done, length, haskey, keys, values, getindex

## Overloads for transparent iteration on `Exponent`, `Polyomial`, `Point`

## Exponent iterator
start(expo::Exponent) = start(expo.expo)
next(expo::Exponent, state) = next(expo.expo, state)
done(expo::Exponent, state) = done(expo.expo, state)
length(expo::Exponent) = length(expo.expo)
haskey(expo::Exponent, key) = haskey(expo.expo, key)
keys(expo::Exponent) = keys(expo.expo)
values(expo::Exponent) = values(expo.expo)
getindex(expo::Exponent, var::Variable) = expo.expo[var]


## Polynomial
start(poly::Polynomial) = start(poly.poly)
next(poly::Polynomial, state) = next(poly.poly, state)
done(poly::Polynomial, state) = done(poly.poly, state)
length(poly::Polynomial) = length(poly.poly)
haskey(poly::Polynomial, key) = haskey(poly.poly, key)
keys(poly::Polynomial) = keys(poly.poly)
values(poly::Polynomial) = values(poly.poly)
getindex(poly::Polynomial, expo::Exponent) = poly.poly[expo]


## Point
start(pt::Point) = start(pt.coords)
next(pt::Point, state) = next(pt.coords, state)
done(pt::Point, state) = done(pt.coords, state)
length(pt::Point) = length(pt.coords)
haskey(pt::Point, key) = haskey(pt.coords, key)
keys(pt::Point) = keys(pt.coords)
values(pt::Point) = values(pt.coords)
getindex(pt::Point, var::Variable) = pt.coords[var]
