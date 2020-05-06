@testset "FM90" begin

    x_inv_microns = [3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 9.0, 10.0]
    wave = 1e4 ./ x_inv_microns

    model = FM90()
    model1 = FM90([0.10 0.70 3.23 0.41])

    # Test out of bounds
    bad_waves = [100, 4e4]
    @test model.(bad_waves) == zeros(length(bad_waves))

    ref_values = [2.9829317, 4.1215415, 6.4135842, 5.6574243, 4.7573250, 5.4905843, 9.2853567, 12.462238]
    @test model.(wave) ≈ ref_values rtol = 1e-4
    @test model.(wave) == model1.(wave)

    # uncertainties
    noise = randn(length(wave)) .* 0.01
    wave_unc = wave .± noise
    reddening = @inferred broadcast(model, wave_unc)
    reddening1 = @inferred broadcast(model1, wave_unc)
    @test Measurements.value.(reddening) ≈ ref_values rtol = 1e-4
    @test Measurements.value.(reddening1) ≈ ref_values rtol = 1e-4

    # Unitful
    wave_u = wave * u"angstrom"
    reddening = @inferred broadcast(model, wave_u)
    reddening1 = @inferred broadcast(model1, wave_u)
    @test eltype(reddening) <: Gain
    @test eltype(reddening1) <: Gain
    @test ustrip.(reddening) ≈ ref_values rtol = 1e-4
    @test ustrip.(reddening1) ≈ ref_values rtol = 1e-4
end
