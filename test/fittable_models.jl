@testset "FM90" begin

    x_inv_microns = [3.5 4.0 4.5 5.0 6.0 7.0 9.0 10.0]
    wave = 1e4 ./ x_inv_microns

    model = FM90()

    # Test out of bounds
    bad_waves = [100, 4e4]
    @test model.(bad_waves) == zeros(length(bad_waves))

    ref_values = [2.9829317 4.1215415 6.4135842 5.6574243 4.7573250 5.4905843 9.2853567 12.462238]
    @test model.(wave) ≈ ref_values rtol = 1e-4

    # uncertainties
    noise = randn(length(wave)) .* 0.01
    wave_unc = wave .± noise
    reddening = model.(wave)
    @test Measurements.value.(reddening) ≈ ref_values rtol = 1e-4

    # # Unitful
    wave_u = wave * u"angstrom"
    reddening = model.(wave_u)
    @test eltype(reddening) <: Gain
    @test ustrip.(reddening) ≈ ref_values rtol = 1e-4
end
