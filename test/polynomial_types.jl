@testset "evaluation tests" begin
    x = MPC.Variable("x", Complex)
    y = MPC.Variable("y", Complex)
    z = MPC.Variable("z", Real)
    b = MPC.Variable("b", Bool)

    ## Polynomial algebra
    pt = Point(SortedDict(x=>2, y=>1+im, z=>0+im, b=>2.1))
    # print(pt)

    p1 = Base.power_by_squaring(x, 2) + 3*y + conj(x) + conj(z) + b
    evaluate(p1, pt) == 10 + 3im

    p2 = Base.power_by_squaring(y,6) - Base.power_by_squaring(y,6) + (-y*x*b) / 4 + π*conj(b)


    ## Poly operators
    @test evaluate(x, pt) == 2
    @test evaluate(y, pt) == 1+im
    @test evaluate(conj(y), pt) == 1-im
    @test evaluate(z, pt) == 0
    @test evaluate(b, pt) == 1

    @test evaluate(x*y, pt) == 2+2im
    @test evaluate(p2, pt) == (-1-im)/2 + π

    @test evaluate(real(y), pt) == 1
    @test evaluate(imag(y), pt) == 1

    @test evaluate(x*y*conj(y), pt) == 2*(1+im)*(1-im)
    @test evaluate(x*y + conj(y), pt) == 3+im
    @test evaluate(x*y + 1 + 1im, pt) == 3+3im


    ## Point algebra
    pt = Point(SortedDict(x=>2, y=>1+im, z=>0+im, b=>2.1))
    pt1 = Point(SortedDict(z=>0+im))
    pt2 = Point(SortedDict(x=>3, y=>-1+2im, b=>-1))
    @test merge(pt, pt1, pt2) == Point(SortedDict(x=>5, y=>3im, z=>0+im, b=>1))

    @test pt + pt1 - 3*pt2 == Point(SortedDict(x=>-7, y=>4-5im, z=>0+im, b=>1))

    @test norm(pt2, Inf) == 3
    @test norm(pt2, 1) == 6
    @test norm(pt2, 2) == √(14)

end
