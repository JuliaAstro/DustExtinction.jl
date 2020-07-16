@testset "G16" begin

    function get_fA_1_ref()
        # From Fitzpatrick (1999) Table 3
        x_invum = [0.377, 0.820, 1.667, 1.828, 2.141, 2.433, 3.704, 3.846]
        ref_values = [0.265, 0.829, 2.688, 3.055, 3.806, 4.315, 6.265, 6.591] ./ 3.1
        return (x_invum, ref_values)
    end

    # Test Rv=3.1, fA=1.0 mixture
    x_invum, ref_values = get_fA_1_ref()
    model = G16(Rv=3.1, f_A=1.0)
    wave = 1e4 ./ x_invum
    @test model.(wave) ≈ ref_values rtol = 2e-3

    # Test RV=2.74, f_A=0 mixture
    tmodel = G16(Rv=2.74, f_A=0.0)
    gmodel = DustExtinction.G03_SMCBar()
    x_invum = gmodel.obsdata_x
    wave = 1e4 ./ x_invum
    ref_values = gmodel.obsdata_axav
    tolerance = gmodel.obsdata_tolerance
    @test collect(tmodel.(wave)) ≈ collect(ref_values) rtol = tolerance

    # Test out of bounds
    model = G16()
    bad_waves = [100, 4e4]
    @test model.(bad_waves) == zeros(length(bad_waves))

    # uncertainties
    x_invum, ref_values = get_fA_1_ref()
    wave = 1e4 ./ x_invum
    model = G16(Rv=3.1, f_A=1.0)
    noise = rand(length(wave)) .* 0.01
    wave_unc = wave .± noise
    reddening = map(w -> @uncertain(model(w)), wave_unc)
    @test Measurements.value.(reddening) ≈ ref_values rtol = 1e-3

    # Unitful
    wave_u = wave * u"angstrom"
    reddening = @inferred map(model, wave_u)
    @test eltype(reddening) <: Gain
    @test ustrip.(reddening) ≈ ref_values rtol = 0.016
end
