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

    refcoords = [0.      0. ;
                 π / 2.  0. ;
                 π       0. ;
                 3π / 2. 0. ;
                 0.1     0.1;
                 0.1    -0.1;
    ]

    refebv = [100.0270,
              2.6185,
              1.4182,
              3.4194,
              0.7949,
              0.5680]

    dustmap = SFD98Map()
    for i = 1:size(refcoords, 1)
        l, b = refcoords[i, :]
        @test dustmap(l, b) ≈ refebv[i] rtol = 0.02
    end
    # vectorized
    l, b = eachcol(refcoords)
    @test dustmap.(l, b) ≈ refebv rtol = 0.02

    # deprecations
    @test_deprecated ebv_galactic(dustmap, l[1], b[1])
    @test_deprecated ebv_galactic(dustmap, l, b)

    @testset "Measurements" begin
        refcoords .± 0.01
        output = dustmap.(eachcol(refcoords)...)
        @test Measurements.values.(output) ≈ refebv rtol = 0.02
    end

    @testset "Unitful" begin
        l, b = eachcol(refcoords * u"rad")
        output = dustmap.(l, b)
        @test eltype(output) <: Gain
        @test ustrip.(output) ≈ refebv rtol = 0.02
    end
end
