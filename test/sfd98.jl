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
