@testset "cplx2real - Exponent" begin
    z1 = MPC.Variable("z1", Complex)
    z2 = MPC.Variable("z2", Complex)
    x = MPC.Variable("x", Real)
    b = MPC.Variable("b", Bool)

    expo_cplx = Exponent(SortedDict(z1=>Degree(1,0),
                                    x=>Degree(4,0),
                                    b=>Degree(1,0)))

    z1_re = MPC.Variable("z1_Re", Real)
    z1_im = MPC.Variable("z1_Im", Real)

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
end

@testset "cplx2real - Polynomial" begin
    z1 = MPC.Variable("z1", Complex)
    z2 = MPC.Variable("z2", Complex)
    x = MPC.Variable("x", Real)
    b = MPC.Variable("b", Bool)

    expo_cplx = Exponent(SortedDict(z1=>Degree(1,0),
                                    x=>Degree(4,0),
                                    b=>Degree(1,0)))
    p = Polynomial(SortedDict{Exponent, Number}(expo_cplx=>3+5im,
                                                Exponent()=>5im))

    # Corresponding real quantities
    z1_re = MPC.Variable("z1_Re", Real)
    z1_im = MPC.Variable("z1_Im", Real)
    expo_real = Exponent(SortedDict(z1_re=>Degree(1,0),
                                    x=>Degree(4,0),
                                    b=>Degree(1,0)))
    expo_imag = Exponent(SortedDict(z1_im=>Degree(1,0),
                                    x=>Degree(4,0),
                                    b=>Degree(1,0)))

    preal, pimag = cplx2real(p)
    @test preal == Polynomial(SortedDict{Exponent,Number}(expo_real=>3.0,
                                                          expo_imag=>-5.0))
    @test pimag == Polynomial(SortedDict{Exponent,Number}(Exponent()=>5,
                                                          expo_imag=>3.0,
                                                          expo_real=>5.0))

end

@testset "cplx2real - Point" begin
    z1 = MPC.Variable("z1", Complex)
    z2 = MPC.Variable("z2", Complex)
    x = MPC.Variable("x", Real)
    b = MPC.Variable("b", Bool)

    pt_c = Point(SortedDict(z1=>1+2im,
                            z2=>3im,
                            x=>3.5,
                            b=>1))

    z1_Re = MPC.Variable("z1_Re", Real)
    z1_Im = MPC.Variable("z1_Im", Real)
    z2_Im = MPC.Variable("z2_Im", Real)

    pt_re = Point(SortedDict(z1_Re=>1,
                             z1_Im=>2,
                             z2_Im=>3,
                             x=>3.5,
                             b=>1))

    @test cplx2real(pt_c) == pt_re

    @test real2cplx(pt_re) == pt_c
end
