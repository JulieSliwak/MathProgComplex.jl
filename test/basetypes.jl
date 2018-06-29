@testset "Variable type" begin
    d = Degree(2,3)

    z = Variable("z", Complex)
    x = Variable("x", Real)
    b = Variable("b", Bool)

    @test iscomplex(z) == true
    @test isreal(x) == true
    @test isbool(b) == true

    show(d)
    show(z)
    println()
end


@testset "Exponent type" begin
    z = Variable("z", Complex)
    x = Variable("x", Real)
    
    one = Exponent()
    zexpo = Exponent(z)
    expo = Exponent(SortedDict(z=>Degree(3,2),
                               x=>Degree(2,0)))
    
    @test zexpo == Exponent(SortedDict(z=>Degree(1,0)))

    show(expo)
    println()
end

@testset "Polynomial type" begin
    emptypoly = Polynomial()

    a = Variable("a", Complex)
    b = Variable("b", Real)

    p_obj = abs2(a) + abs2(b) + 2
    p_cstr1 = 3a + b + 2
    p_cstr2 = abs2(b) + 5a*b + 2
    
    io = IOBuffer()
    print(io, p_obj)
    @test String(take!(io)) == "2 + conj(a) * a + b^2"
    print(io, p_cstr1)
    @test String(take!(io)) == "2 + (3.0)*a + b"
    print(io, p_cstr2)
    @test String(take!(io)) == "2 + (5.0)*a * b + b^2"


end