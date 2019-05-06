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

    for i = 1:length(wave)
        @test ccm89(wave[i], 3.1) ≈ ref_values[i] rtol = 0.016
    end

    # Now test for dot operation
    @test ccm89.(wave, 3.1) ≈ ref_values rtol=0.016

    # Test errors
    bad_waves = [
        100.:13000.,
        33400.:40000.
    ]
    for bad_wave in bad_waves
        @test_throws ErrorException ccm89.(bad_wave, 3.1)
    end
end