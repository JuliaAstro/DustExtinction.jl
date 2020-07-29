using Measurements
using DustExtinction: aa_to_invum,
                      g03_invum,
                      g16_invum,
                      G03_SMCBar,
                      G16

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
    model = G16(Rv=3.1, f_A=1.0)
    bad_waves = [100, 4e4]
    @test model.(bad_waves) == zeros(length(bad_waves))
    @test_throws ErrorException g16_invum(aa_to_invum(bad_waves[1]), 3.1, 1.0)
    @test_throws ErrorException g16_invum(aa_to_invum(bad_waves[2]), 3.1, 1.0)

    # uncertainties
    x_invum, ref_values = get_fA_1_ref()
    wave = 1e4 ./ x_invum
    model = G16(Rv=3.1, f_A=1.0)
    noise = rand(length(wave)) .* 0.01
    wave_unc = wave .± noise
    reddening = map(w -> @uncertain(model(w)), wave_unc)
    @test Measurements.value.(reddening) ≈ ref_values rtol = 1e-3

    # Unitful
    model = G16(Rv=3.1, f_A=1.0)
    x_invum, ref_values = get_fA_1_ref()
    wave = 1e4 ./ x_invum
    wave_u = wave * u"angstrom"
    reddening = @inferred map(model, wave_u)
    @test eltype(reddening) <: Gain
    @test ustrip.(reddening) ≈ ref_values rtol = 0.016

end

@testset "G03_SMCBar" begin

    # Test output
    model = G03_SMCBar()
    x, y, tolerance = model.obsdata_x, model.obsdata_axav, model.obsdata_tolerance
    w = 1e4 ./ x
    @test model.(collect(w)) ≈ collect(y) rtol = tolerance

    # Test out of bounds
    model = G03_SMCBar(Rv=2.74)
    bad_waves = [100, 4e4]
    @test model.(bad_waves) == zeros(length(bad_waves))
    @test_throws ErrorException g03_invum(aa_to_invum(bad_waves[1]), 3.1)
    @test_throws ErrorException g03_invum(aa_to_invum(bad_waves[2]), 3.1)

    # uncertainties
    model = G03_SMCBar(Rv=2.74)
    x, y, tolerance = model.obsdata_x, model.obsdata_axav, model.obsdata_tolerance
    wave = 1e4 ./ x
    noise = rand(length(wave)) .* 0.01
    wave_unc = wave .± noise
    reddening = map(w -> @uncertain(model(w)), wave_unc)
    @test Measurements.value.(reddening) ≈ collect(y) rtol = tolerance

    # Unitful
    model = G03_SMCBar(Rv=2.74)
    x, y, tolerance = model.obsdata_x, model.obsdata_axav, model.obsdata_tolerance
    wave = 1e4 ./ x
    wave_u = collect(wave) * u"angstrom"
    reddening = @inferred map(model, wave_u)
    @test eltype(reddening) <: Gain
    @test ustrip.(reddening) ≈ collect(y) rtol = tolerance
end
