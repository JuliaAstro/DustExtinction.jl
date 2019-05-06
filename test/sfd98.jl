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

    for i = 1:length(wave)
        @test od94(wave[i], 3.1) ≈ ref_values[i] rtol = 0.0051
    end

    # Now test for dot operation
    @test od94.(wave, 3.1) ≈ ref_values rtol=0.0051

    # Test errors
    bad_waves = [
        100.:13000.,
        33400.:40000.
    ]
    for bad_wave in bad_waves
        @test_throws ErrorException od94.(bad_wave, 3.1)
    end
end

@testset "SFD98 dust maps" begin

    # refebv obtained using http://irsa.ipac.caltech.edu/applications/DUST/
    # and manually inserting the following lines and reading off the values for
    # SFD (1998) reference pixel:
    #
    # 0. 0. gal
    # 90. 0. gal
    # 180. 0. gal
    # 270. 0. gal
    # 5.729577951308233 5.729577951308233 gal
    # 5.729577951308233 -5.729577951308233 gal
    #
    # Note that values are not expected to agree exactly because (1) IRSA
    # seems to not do linear interpolation and (2) IRSA seems to use a
    # resampled map (pixel values reported by IRSA don't seem to match any
    # pixel values in original maps).

    refcoords = [(0., 0.),
                (pi / 2., 0.),
                (pi, 0.),
                (3pi / 2., 0.),
                (0.1, 0.1),
                (0.1, -0.1)]

    refebv = [100.0270,
            2.6185,
            1.4182,
            3.4194,
            0.7949,
            0.5680]

    if haskey(ENV, "SFD98_DIR")
        dustmap = SFD98Map()
        for i = 1:length(refcoords)
            l, b = refcoords[i]
            @test ebv_galactic(dustmap, l, b) ≈ refebv[i] rtol = 0.02
        end

    else
        println("Skipping SFD98Map test because \$SFD98_DIR not defined.")
    end
end