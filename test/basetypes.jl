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
    @test Exponent(SortedDict(z=>Degree(3,2),
                              x=>Degree(0,0))) == Exponent(SortedDict(z=>Degree(3,2)))

    show(one); println()
    show(expo); println()
end

@testset "Polynomial type" begin
    emptypoly = Polynomial()

    a = Variable("a", Complex)
    b = Variable("b", Real)

    one = Exponent()
    expoa2 = Exponent(SortedDict(a=>Degree(1,1)))
    expob2 = Exponent(SortedDict(b=>Degree(1,0)))
    expoa = Exponent(SortedDict(a=>Degree(1,0)))
    expob = Exponent(SortedDict(b=>Degree(1,0)))

    p_obj = Polynomial(SortedDict{Exponent, Number}(expoa2=>1.0,
                                                    expob2=>3.5,
                                                    expoa=>5+3im,
                                                    expob=>0,
                                                    one=>-1))
    
    io = IOBuffer()
    print(io, p_obj)
    @test String(take!(io)) == "-1 + (5 + 3im)*a + conj(a) * a"
end