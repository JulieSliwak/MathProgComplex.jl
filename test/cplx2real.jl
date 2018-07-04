@testset "cplx2real - Exponent" begin
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

    preal, pimag = cplx2real(expo_cplx)
    @test preal == Polynomial(SortedDict{Exponent, Number}(expo_real=>1.0))
    @test pimag == Polynomial(SortedDict{Exponent, Number}(expo_imag=>1.0))

end
