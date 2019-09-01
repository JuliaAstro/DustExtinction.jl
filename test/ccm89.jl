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
    x_inv_microns = [2.78; 2.27; 1.82; 1.43; 1.11; 0.80; 0.63; 0.46]
    wave = 1.e4 ./ x_inv_microns

    # A(lambda)/A(V) for R_V = 3.1 from Table 3 of CCM '89
    ref_values = [1.569, 1.337, 1.000, 0.751, 0.479, 0.282, 0.190, 0.114]

    # Now test for dot operation
    @test ccm89.(wave, 3.1) ≈ ref_values rtol = 0.016

    # Test deprecated array syntax
    @test_deprecated ccm89(wave)

    # Test errors
    bad_waves = [
        100.:13000.,
        33400.:40000.
    ]
    for bad_wave in bad_waves
        @test_throws ErrorException ccm89.(bad_wave, 3.1)
    end

    @testset "uncertainties" begin        
        noise = randn(length(wave)) .* 0.01
        wave_unc = wave .± noise
        reddening = ccm89.(wave_unc, 3.1)
        @test Measurements.value.(reddening) ≈ ref_values rtol = 0.016
    end 

    @testset "unitful" begin
        wave_u = wave * u"angstrom"
        reddening = ccm89.(wave_u, 3.1)
        @test eltype(reddening) <: Gain
        @test ustrip.(reddening) ≈ ref_values rtol = 0.016
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
            
                # Now test for dot operation
    @test od94.(wave, 3.1) ≈ ref_values rtol = 0.016
            
                # Test deprecated array syntax
    @test_deprecated od94(wave)
            
                # Test errors
    bad_waves = [
                    100.:13000.,
                    33400.:40000.
                ]
    for bad_wave in bad_waves
        @test_throws ErrorException od94.(bad_wave, 3.1)
    end

    @testset "uncertainties" begin        
        noise = randn(length(wave)) .* 10
        wave_unc = wave .± noise
        reddening = od94.(wave_unc, 3.1)
        @test Measurements.value.(reddening) ≈ ref_values rtol = 0.016
    end
    @testset "unitful" begin
        wave_u = wave * u"angstrom"
        reddening = od94.(wave_u, 3.1)
        @test eltype(reddening) <: Gain
        @test ustrip.(reddening) ≈ ref_values rtol = 0.016
    end
end
