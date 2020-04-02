using Measurements

@testset "ccm89" begin

    # NOTE: Test is only to precision of 0.015 because there is a discrepancy
    # of 0.014 for the B band wavelength of unknown origin (and up to 0.002 in
    # other bands).
    #
    # Note that a and b can be obtained with:
    # b = ccm89(wave, 0.)
    # a = ccm89(wave, 1.) - b
    #
    # These differ from the values tablulated in the original paper.
    # Could be due to floating point errors in the original paper?

    # U, B, V, R, I, J, H, K band effective wavelengths from CCM '89 table 3
    x_inv_microns = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.6, 4.0, 3.0, 2.78, 2.27, 1.82, 1.43, 1.11, 0.8, 0.63, 0.46]
    wave = 1e4 ./ x_inv_microns

    # A(lambda)/A(V) for different R_V from Table 3 of CCM '89
    ref_values = Dict(3.1 => [5.23835484, 4.13406452, 3.33685933, 2.77962453, 2.52195399, 2.84252644, 3.18598916, 2.31531711, 1.64254927, 1.56880904, 1.32257836, 1.0, 0.75125994, 0.4780346, 0.28206957, 0.19200814, 0.11572348],
        2.0 => [9.407, 7.3065, 5.76223881, 4.60825807, 4.01559036, 4.43845534, 4.93952892, 3.39275574, 2.068771, 1.9075018, 1.49999733, 1.0, 0.68650255, 0.36750326, 0.21678862, 0.14757062, 0.08894094],
        3.0 => [5.491, 4.32633333, 3.48385202, 2.8904508, 2.6124774, 2.9392494, 3.2922643, 2.38061642, 1.66838089, 1.58933588, 1.33333103, 1.0, 0.74733525, 0.47133573, 0.27811315, 0.18931496, 0.11410029],
        4.0 => [3.533, 2.83625, 2.34465863, 2.03154717, 1.91092092, 2.18964643, 2.46863199, 1.87454675, 1.46818583, 1.43025292, 1.24999788, 1.0, 0.7777516, 0.52325196, 0.30877542, 0.21018713, 0.12667997],
        5.0 => [2.3582, 1.9422, 1.66114259, 1.51620499, 1.48998704, 1.73988465, 1.97445261, 1.57090496, 1.3480688, 1.33480314, 1.19999799, 1.0, 0.79600141, 0.5544017, 0.32717278, 0.22271044, 0.13422778],
        6.0 => [1.575, 1.34616667, 1.20546523, 1.17264354, 1.20936444, 1.44004346, 1.64499968, 1.36847709, 1.26799077, 1.27116996, 1.16666472, 1.0, 0.80816794, 0.5751682, 0.33943769, 0.23105931, 0.13925965])

    # Test defaults
    @test @inferred(broadcast(ccm89, wave)) ≈ ref_values[3.1] rtol = 0.016
    @test_deprecated ccm89(wave)

    for rv in [2.0, 3.0, 3.1, 4.0, 5.0, 6.0]
        output = @inferred broadcast(ccm89, wave, rv)
        @test output ≈ ref_values[rv] rtol = 0.016

        # Test deprecated array syntax
        @test_deprecated ccm89(wave, rv)

        bad_waves = [100, 4e4]
        @test @inferred(broadcast(ccm89, bad_waves, rv)) == zeros(length(bad_waves))

        # uncertainties
        noise = randn(length(wave)) .* 0.01
        wave_unc = wave .± noise
        reddening = @inferred broadcast(ccm89, wave_unc, rv)
        @test Measurements.value.(reddening) ≈ ref_values[rv] rtol = 0.016

        # Unitful
        wave_u = wave * u"angstrom"
        reddening = @inferred broadcast(ccm89, wave_u, rv)
        @test eltype(reddening) <: Gain
        @test ustrip.(reddening) ≈ ref_values[rv] rtol = 0.016
    end
end

@testset "od94" begin
    # NOTE: The tabulated values go to 0.001, but the test is only for matching
    # at the 0.005 level, because there is currently a discrepancy up to 0.0047
    # of unknown origin.

    # Tests od94() at Rv = 3.1 against the widely used values tabulated in
    # Schlegel, Finkbeiner and Davis (1998)
    # http://adsabs.harvard.edu/abs/1998ApJ...500..525S

    # This is tested by evaluating the extinction curve at a (given)
    # effective wavelength, since these effective wavelengths:
    # "... represent(s) that wavelength on the extinction curve
    # with the same extinction as the full passband."

    # The test does not include UKIRT L' (which, at 3.8 microns) is
    # beyond the range of wavelengths allowed by the function
    # or the APM b_J filter which is defined in a non-standard way.

    # The SFD98 tabulated values go to 1e-3, so we should be able to match at
    # that level.

    wave = [3372., 4404., 5428., 6509., 8090.,
            3683., 4393., 5519., 6602., 8046.,
            12660., 16732., 22152.,
            5244., 6707., 7985., 9055.,
            6993.,
            3502., 4676., 4127.,
            4861., 5479.,
            3546., 4925., 6335., 7799., 9294.,
            3047., 4711., 5498.,
            6042., 7068., 8066.,
            4814., 6571., 8183.]

    ref_values = [1.664, 1.321, 1.015, 0.819, 0.594,
                1.521, 1.324, 0.992, 0.807, 0.601,
                0.276, 0.176, 0.112,
                1.065, 0.793, 0.610, 0.472,
                0.755,
                1.602, 1.240, 1.394,
                1.182, 1.004,
                1.579, 1.161, 0.843, 0.639, 0.453,
                1.791, 1.229, 0.996,
                0.885, 0.746, 0.597,
                1.197, 0.811, 0.580]

    reddening = @inferred broadcast(od94, wave, 3.1)
    @test reddening ≈ ref_values rtol = 0.016

    # Test deprecated array syntax
    @test_deprecated od94(wave)

    # Test out of bounds
    bad_waves = [100, 4e4]
    @test @inferred(broadcast(od94, bad_waves, 3.1)) == zeros(length(bad_waves))

    @testset "uncertainties" begin
        noise = randn(length(wave)) .* 10
        wave_unc = wave .± noise
        reddening = @inferred broadcast(od94, wave_unc, 3.1)
        @test Measurements.value.(reddening) ≈ ref_values rtol = 0.016
    end
    @testset "unitful" begin
        wave_u = wave * u"angstrom"
        reddening = @inferred broadcast(od94, wave_u, 3.1)
        @test eltype(reddening) <: Gain
        @test ustrip.(reddening) ≈ ref_values rtol = 0.016
    end
end

@testset "cal00" begin

    refwave = [3090.90909091,  4561.61616162,  6872.72727273,  9604.04040404,
    14646.46464646, 14646.46464646, 15486.86868687, 18218.18181818,
    20529.29292929, 20949.49494949]

    refmag = Dict(3.1 => [ 1.88010678,  1.27137591,  0.70513192,  0.33600273,  0.01622915,
    0.01622915, -0.01682162, -0.10317769, -0.15830055, -0.16701622],
        2.4 => [ 2.13680458,  1.35052722,  0.61912873,  0.14233687, -0.27070401,
    -0.27070401, -0.3133946 , -0.42493784, -0.49613821, -0.50739595])

    # test defaults
    @test cal00.(refwave, 3.1) ≈ refmag[3.1]

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

@testset "vcg04" begin

    x_inv_microns = [8.0, 7.0, 6.0, 5.0, 4.6, 4.0, 3.4]
    wave = 1e4 ./ x_inv_microns

    ref_values = Dict(3.1 => [3.36528, 2.84166, 2.58283, 2.88248, 3.25880, 2.43315, 2.00025],
        2.0 => [5.20767, 4.25652, 3.74640, 4.16150, 4.73050, 3.33399, 2.54668],
        3.0 => [3.47694, 2.92741, 2.65335, 2.96000, 3.34799, 2.48775, 2.03337],
        4.0 => [2.61157, 2.26285, 2.10683, 2.35925, 2.65674, 2.06463, 1.77671],
        5.0 => [2.09235, 1.86411, 1.77892, 1.99880, 2.24199, 1.81076, 1.622711],
        6.0 => [1.74620, 1.59829, 1.56031, 1.75850, 1.96549, 1.64151, 1.52005])

    @test @inferred(broadcast(vcg04, wave)) ≈ ref_values[3.1] rtol = 0.016
    @test_deprecated vcg04(wave)

    for rv in [2.0, 3.0, 3.1, 4.0, 5.0, 6.0]
        output = @inferred broadcast(vcg04, wave, rv)
        @test output ≈ ref_values[rv] rtol = 0.016

        # Test deprecated array syntax
        @test_deprecated vcg04(wave, rv)

        bad_waves = [100, 4e4]
        @test @inferred(broadcast(vcg04, bad_waves, rv)) == zeros(length(bad_waves))

        # uncertainties
        noise = randn(length(wave)) .* 0.01
        wave_unc = wave .± noise
        reddening = @inferred broadcast(vcg04, wave_unc, rv)
        @test Measurements.value.(reddening) ≈ ref_values[rv] rtol = 0.016

        # Unitful
        wave_u = wave * u"angstrom"
        reddening = @inferred broadcast(vcg04, wave_u, rv)
        @test eltype(reddening) <: Gain
        @test ustrip.(reddening) ≈ ref_values[rv] rtol = 0.016
    end
end

@testset "gcc09" begin

    x_inv_microns = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.6, 4.0, 3.4]
    wave = 1e4 ./ x_inv_microns

    ref_values = Dict(3.1 => [5.23161, 4.20810, 3.45123, 2.92264, 2.61283, 2.85130, 3.19451, 2.34301, 1.89256],
        2.0 => [10.5150, 8.07274, 6.26711, 5.00591, 4.24237, 4.42844, 4.99482, 3.42585, 2.59322],
        3.0 => [5.55181, 4.44232, 3.62189, 3.04890, 2.71159, 2.94688, 3.30362, 2.40863, 1.93502],
        4.0 => [3.07020, 2.62711, 2.29927, 2.07040, 1.94621, 2.20610, 2.45801, 1.90003, 1.60592],
        5.0 => [1.58123, 1.53798, 1.50571, 1.48330, 1.48697, 1.76164, 1.95065, 1.59486, 1.40846],
        6.0 => [0.588581, 0.811898, 0.976660, 1.09190, 1.18082, 1.46533, 1.61241, 1.39142, 1.27682])


    @test @inferred(broadcast(gcc09, wave)) ≈ ref_values[3.1] rtol = 0.016
    @test_deprecated gcc09(wave)

    for rv in [2.0, 3.0, 3.1, 4.0, 5.0, 6.0]
        output = @inferred broadcast(gcc09, wave, rv)
        @test output ≈ ref_values[rv] rtol = 0.016

        # Test deprecated array syntax
        @test_deprecated gcc09(wave, rv)

        bad_waves = [100, 4e4]
        @test @inferred(broadcast(gcc09, bad_waves, rv)) == zeros(length(bad_waves))

        # uncertainties
        noise = randn(length(wave)) .* 0.01
        wave_unc = wave .± noise
        reddening = @inferred broadcast(gcc09, wave_unc, rv)
        @test Measurements.value.(reddening) ≈ ref_values[rv] rtol = 0.016

        # Unitful
        wave_u = wave * u"angstrom"
        reddening = @inferred broadcast(gcc09, wave_u, rv)
        @test eltype(reddening) <: Gain
        @test ustrip.(reddening) ≈ ref_values[rv] rtol = 0.016
    end
end
