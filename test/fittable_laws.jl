@testset "FM90" begin

    x_inv_microns = [3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 9.0, 10.0]
    wave = 1e4 ./ x_inv_microns

    # Construction
    model = FM90()
    model1 = FM90([0.10 0.70 3.23 0.41])
    @test_throws ErrorException("`x0` must be ≥ 0, got -2.0") FM90(x0=-2.0)
    @test_throws ErrorException("`gamma` must be ≥ 0, got -2.0") FM90(gamma=-2.0)

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

@testset "P92" begin

    x_inv_microns = [0.21, 0.29, 0.45, 0.61, 0.80, 1.11, 1.43, 1.82, 2.27, 2.50, 2.91, 3.65, 4.00, 4.17, 4.35,
                     4.57, 4.76, 5.00, 5.26, 5.56, 5.88, 6.25, 6.71, 7.18, 7.60, 8.00, 8.50, 9.00, 9.50, 10.00]

    wave = 1e4 ./ x_inv_microns

    MW_exvebv = [-3.02, -2.91, -2.76, -2.58, -2.23, -1.60, -0.78, 0.00, 1.00, 1.30, 1.80, 3.10, 4.19, 4.90, 5.77,
                  6.57, 6.23, 5.52, 4.90, 4.65, 4.60, 4.73, 4.99, 5.36, 5.91, 6.55, 7.45, 8.45, 9.80, 11.30]

    Rv = 3.08
    ref_values = MW_exvebv ./ Rv .+ 1

    model = P92()
    model1 = P92(FUV_b = 4)

    # Test out of bounds
    bad_waves = [9, 1e8]
    @test model.(bad_waves) == zeros(length(bad_waves))
    @test model1.(bad_waves) == zeros(length(bad_waves))

    # testing main part
    @test model.(wave) ≈ ref_values rtol = 0.25 atol = 0.01
    @test model1.(wave) ≈ ref_values rtol = 0.25 atol = 0.01

    # uncertainties
    noise = randn(length(wave)) .* 0.01
    wave_unc = wave .± noise
    reddening = @inferred broadcast(model, wave_unc)
    reddening1 = @inferred broadcast(model1, wave_unc)
    @test Measurements.value.(reddening) ≈ ref_values rtol = 0.25 atol = 0.01
    @test Measurements.value.(reddening1) ≈ ref_values rtol = 0.25 atol = 0.01

    # Unitful
    wave_u = wave * u"angstrom"
    reddening = @inferred broadcast(model, wave_u)
    reddening1 = @inferred broadcast(model1, wave_u)
    @test eltype(reddening) <: Gain
    @test eltype(reddening1) <: Gain
    @test ustrip.(reddening) ≈ ref_values rtol = 0.25 atol = 0.01
    @test ustrip.(reddening1) ≈ ref_values rtol = 0.25 atol = 0.01
end
