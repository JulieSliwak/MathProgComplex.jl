using NamedArrays, BenchmarkTools

I = [randstring(4) for i=1:100]
J = [randstring(4) for j=1:100]

d1 = Dict((i,j) => rand() for i in I, j in J)
d2 = Dict((i,j) => rand() for i in I, j in J)

d3 = Dict((i,j) => d1[i,j]*d2[i,j] for i in I, j in J)
d4 = Dict{Tuple{String,String}, Float64}()
for i in I, j in J
    d4[i,j] = d1[i,j]*d2[i,j]
end

println("\nDict comprehension\n", @benchmark d3 = Dict((i,j) => d1[i,j]*d2[i,j] for i in I, j in J))
println("\nDict assignment\n", @benchmark begin
    d4 = Dict{Tuple{String,String}, Float64}()
    for i in I, j in J
        d4[i,j] = d1[i,j]*d2[i,j]
    end
end)

n1 = NamedArray([rand() for i in I, j in J], (I,J))
n2 = NamedArray([rand() for i in I, j in J], (I,J))

n3 = NamedArray([n1[i,j]*n2[i,j] for i in I, j in J], (I,J))
n4 = NamedArray(zeros(100,100), (I,J))
for i in I, j in J
    n4[i,j] = n1[i,j]*n2[i,j]
end
println("\nNamedArray comprehension\n", @benchmark n3 = NamedArray([n1[i,j]*n2[i,j] for i in I, j in J], (I,J)))
println("\nNamedArray assignment\n", @benchmark begin
    n4 = NamedArray(zeros(100,100), (I,J))
    for i in I, j in J
        n4[i,j] = n1[i,j]*n2[i,j]
    end
end)

b1 = [rand() for i=1:100, j=1:100]
b2 = [rand() for i=1:100, j=1:100]

b3 = [b1[i,j]*b2[i,j] for i=1:100, j=1:100]
println("\nBasic Array comprehension (no name lookup)\n", @benchmark b3 = [b1[i,j]*b2[i,j] for i=1:100, j=1:100])

i, j = I[4], J[4]
println("\nBasic Dict comprehension (no name lookup)\n", @benchmark d3[i, j])
println("\nBasic Dict comprehension (no name lookup)\n", @benchmark d4[i, j])

println("\nNamedArray assignment\n", @benchmark n3[i, j])
println("\nNamedArray assignment\n", @benchmark n4[i, j])

println("\nBasic Array comprehension (no name lookup)\n", @benchmark b3[4, 4])