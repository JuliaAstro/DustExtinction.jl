@testset "cal00" begin

    refwave = [3090.90909091,  4561.61616162,  6872.72727273,  9604.04040404,
    14646.46464646, 14646.46464646, 15486.86868687, 18218.18181818,
    20529.29292929, 20949.49494949]

    refmag = Dict(3.1 => [ 1.88010678,  1.27137591,  0.70513192,  0.33600273,  0.01622915,
    0.01622915, -0.01682162, -0.10317769, -0.15830055, -0.16701622],
        2.4 => [ 2.13680458,  1.35052722,  0.61912873,  0.14233687, -0.27070401,
    -0.27070401, -0.3133946 , -0.42493784, -0.49613821, -0.50739595])

    # test defaults
    @test cal00.(refwave) ≈ refmag[3.1]

    for rv in [2.4, 3.1]
        reddening = @inferred broadcast(cal00, refwave, rv)
        @test reddening ≈ refmag[rv]

        bad_waves = [
            1e2,
            3e4,
            4e4
        ]
        @test @inferred(broadcast(cal00, bad_waves, rv)) == zeros(length(bad_waves))

        # Uncertainties
        noise = randn(length(refwave)) .* 10
        wave_unc = refwave .± noise
        reddening = @inferred broadcast(cal00, refwave, rv)
        @test Measurements.value.(reddening) ≈ refmag[rv]


        # Unitful
        wave_u = refwave * u"angstrom"
        reddening = @inferred broadcast(cal00, wave_u, rv)
        @test eltype(reddening) <: Gain
        @test ustrip.(reddening) ≈ refmag[rv]  
    end
end 
