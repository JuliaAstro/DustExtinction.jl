
@testset "cal00" begin

    refwave = [3090.90909091,  4561.61616162,  6872.72727273,  9604.04040404,
    14646.46464646, 14646.46464646, 15486.86868687, 18218.18181818,
    20529.29292929, 20949.49494949]

    refmag_31 = [ 1.88010678,  1.27137591,  0.70513192,  0.33600273,  0.01622915,
    0.01622915, -0.01682162, -0.10317769, -0.15830055, -0.16701622]

    refmag_24 = [ 2.13680458,  1.35052722,  0.61912873,  0.14233687, -0.27070401,
    -0.27070401, -0.3133946 , -0.42493784, -0.49613821, -0.50739595]

    @test cal00.(refwave) ≈ refmag_31
    @test cal00.(refwave, 2.4) ≈ refmag_24


    bad_waves = [
        0:1200,
        22000:40000
    ]
    for bad_wave in bad_waves
        @test_throws ErrorException cal00.(bad_wave)
    end
end