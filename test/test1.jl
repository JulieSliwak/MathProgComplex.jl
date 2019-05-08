include("..\\src\\MathProgComplex.jl")
using Test

d = Degree(2,3)

z = Variable("z", Complex)
x = Variable("x", Real)
b = Variable("b", Bool)

@test iscomplex(z) == true
@test isreal(x) == true
@test isbool(b) == true

@test hash(Degree(1,3)) == hash(Degree(1,3))
@test hash(Variable("a",Complex)) == hash(Variable("a",Complex))

z = Variable("z", Complex)
x = Variable("x", Real)

one = Exponent()
zexpo = Exponent(z)
expo = Exponent(SortedDict(z=>Degree(3,2),
                           x=>Degree(2,0)))

@test zexpo == Exponent(SortedDict(z=>Degree(1,0)))
@test Exponent(SortedDict(z=>Degree(3,2),
                          x=>Degree(0,0))) == Exponent(SortedDict(z=>Degree(3,2)))


emptypoly = Polynomial()

a = Variable("a", Complex)
b = Variable("b", Real)

one = Exponent()
expoa2 = Exponent(SortedDict(a=>Degree(1,1)))
expob2 = Exponent(SortedDict(b=>Degree(1,0)))
expoa = Exponent(SortedDict(a=>Degree(1,0)))
expob = Exponent(SortedDict(b=>Degree(1,0)))

# println(isless(expoa, expoa2))
# println(isequal(expoa2, expoa2))
# println(isless(expoa2, expoa2))
# println(isless_degree(expoa, expoa2))
# println(isless_degree(expoa2, expoa2))

p_obj = Polynomial(SortedDict{Exponent, Number}(expoa2=>1.0,
                                              expob2=>3.5,
                                              expoa=>5+3im,
                                              expob=>0,
                                              one=>-1))

io = IOBuffer()
print(io, p_obj)
@test String(take!(io)) == "-1 + (5 + 3im)*a + conj(a) * a"


    pt = Point()

    z = Variable("z", Complex)
    x1 = Variable("x1", Real)
    x2 = Variable("x2", Real)
    b = Variable("b", Bool)
    b1 = Variable("b1", Bool)

    pointdict = SortedDict(z =>1+2im,
                          x1 => 3.5,
                          x2 => 0,
                          b => 1,
                          b1 => 3.5)


    pt = Point(pointdict)
    @test pt == Point([z x1 x2 b b1], [1+2im 3.5 0 1 3.5])
    @test length(pt) == 4

    io = IOBuffer()
    print(io, pt)
    @test String(take!(io)) == "b  1\nb1 1\nx1 3.5\nz  1 + 2im\n"

    a = Variable("a", Complex)
    b = Variable("b", Real)

    one = Exponent()
    expoa2 = Exponent(SortedDict(a=>Degree(1,1)))
    expob2 = Exponent(SortedDict(b=>Degree(1,0)))
    expoa = Exponent(SortedDict(a=>Degree(1,0)))
    expob = Exponent(SortedDict(b=>Degree(1,0)))

    p1 = Polynomial(SortedDict{Exponent, Number}(expoa2=>1.0,
                                                 expob2=>3.5,
                                                 expoa=>5+3im,
                                                 expob=>0,
                                                 one=>-1))

    ub_c = 3+5im
    lb_c = -3-2im

    C1 = p1 << ub_c
    C2 = ub_c >> p1
    @test C1 == Constraint(p1, -Inf-Inf*im, ub_c)
    @test C2 == Constraint(p1, -Inf-Inf*im, ub_c)

    @test (lb_c << C1) == Constraint(p1, lb_c, ub_c)
    @test (C2 >> lb_c) == Constraint(p1, lb_c, ub_c)

    C1 = lb_c << p1
    C2 = p1 >> lb_c
    @test C1 == Constraint(p1, lb_c, +Inf+Inf*im)
    @test C2 == Constraint(p1, lb_c, +Inf+Inf*im)

    @test (C1 << ub_c) == Constraint(p1, lb_c, ub_c)
    @test (ub_c >> C2) == Constraint(p1, lb_c, ub_c)

    @test (p1 == ub_c) == Constraint(p1, ub_c, ub_c)
    @test (p1 == ub_c) != (p1 << ub_c)

    io = IOBuffer()
    print(io, Constraint(p1, lb_c, ub_c))
    @test String(take!(io)) == "-3 - 2im < -1 + (5 + 3im)*a + conj(a) * a < 3 + 5im"

    print(io, p1 << ub_c)
    @test String(take!(io)) == "-1 + (5 + 3im)*a + conj(a) * a < 3 + 5im"
    print(io, lb_c << p1)
    @test String(take!(io)) == "-3 - 2im < -1 + (5 + 3im)*a + conj(a) * a"

    print(io, Constraint(p1, ub_c, ub_c))
    @test String(take!(io)) == "-1 + (5 + 3im)*a + conj(a) * a = 3 + 5im"


    a = Variable("a", Complex)
    b = Variable("b", Real)

    one = Exponent()
    expoa2 = Exponent(SortedDict(a=>Degree(1,1)))
    expob2 = Exponent(SortedDict(b=>Degree(1,0)))
    expoa = Exponent(SortedDict(a=>Degree(1,0)))
    expob = Exponent(SortedDict(b=>Degree(1,0)))

    p1 = Polynomial(SortedDict{Exponent, Number}(expoa2=>1.0,
                                                expob2=>3.5,
                                                expoa=>5+3im,
                                                expob=>0,
                                                one=>-1))

    vars = SortedDict{String, Type}("a" => Complex, "b" => Real)
    ctrs = SortedDict{String, Constraint}("ctr1" => ((-3+5im) << p1),
                                          "ctr2" => (p1 == 2.5))
    pb = Problem(p1, ctrs, vars)
    Problem()

    io = IOBuffer()
    print(io, pb)
    @test String(take!(io)) == "▶ variables: a b \n▶ objective: -1 + (5 + 3im)*a + conj(a) * a\n▶ constraints: \n →       ctr1: -3 + 5im < -1 + (5 + 3im)*a + conj(a) * a\n →       ctr2: -1 + (5 + 3im)*a + conj(a) * a = 2.5\n"


    z1 = Variable("z1", Complex)
    z2 = Variable("z2", Complex)
    x = Variable("x", Real)
    b = Variable("b", Bool)

    expo_cplx = Exponent(SortedDict(z1=>Degree(1,0),
                                    x=>Degree(4,0),
                                    b=>Degree(1,0)))

    z1_re = Variable("z1_Re", Real)
    z1_im = Variable("z1_Im", Real)

    expo_real = Exponent(SortedDict(z1_re=>Degree(1,0),
                                    x=>Degree(4,0),
                                    b=>Degree(1,0)))
    expo_imag = Exponent(SortedDict(z1_im=>Degree(1,0),
                                    x=>Degree(4,0),
                                    b=>Degree(1,0)))

    preal, pimag = cplx2real(expo_cplx, 1.)
    @test preal == Polynomial(SortedDict{Exponent, Number}(expo_real=>1.0))
    @test pimag == Polynomial(SortedDict{Exponent, Number}(expo_imag=>1.0))


    expo_cplx = Exponent(SortedDict(z1=>Degree(0,1),
                                    x=>Degree(4,0),
                                    b=>Degree(1,0)))

    preal, pimag = cplx2real(expo_cplx, 1.)
    @test preal == Polynomial(SortedDict{Exponent, Number}(expo_real=>1.0))
    @test pimag == Polynomial(SortedDict{Exponent, Number}(expo_imag=>-1.0))
