using DataDeps
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

    ref_l = [0., π / 2., π, 3π / 2., 0.1, 0.1 ]
    ref_b = [0., 0., 0., 0., 0.1, -0.1]

    refebv = [100.0270,
              2.6185,
              1.4182,
              3.4194,
              0.7949,
              0.5680]

    dustmap = SFD98Map()
    for i = 1:length(refebv)
        @test dustmap(ref_l[i], ref_b[i]) ≈ refebv[i] rtol = 0.02
    end
    # vectorized
    @test dustmap.(ref_l, ref_b) ≈ refebv rtol = 0.02

    # deprecations
    @test_deprecated ebv_galactic(dustmap, ref_l[1], ref_b[1])
    @test_deprecated ebv_galactic(dustmap, ref_l, ref_b)

    @testset "Measurements" begin
        ref_l_m = ref_l .± 0.03
        ref_b_m = ref_b .± 0.01
        output = dustmap.(ref_l_m, ref_b_m)
        @test Measurements.values.(output) ≈ refebv rtol = 0.02
    end

    @testset "SkyCoords" begin
        for i in eachindex(refebv)
            @test dustmap(GalCoords(ref_l[i], ref_b[i])) ≈ refebv[i] rtol = 0.02
            @test dustmap(convert(ICRSCoords, GalCoords(ref_l[i], ref_b[i]))) ≈ refebv[i] rtol = 0.02
        end
        @test dustmap.(convert.(Ref(ICRSCoords), GalCoords.(ref_l, ref_b))) ≈ refebv rtol = 0.02
    end

    @testset "Unitful" begin
        ref_l_u = ref_l * u"rad"
        ref_b_u = ref_b * u"rad"
        output = dustmap.(ref_l_u, ref_b_u)
        @test eltype(output) <: Gain
        @test ustrip.(output) ≈ refebv rtol = 0.02
    end
    mapdir = datadep"sfd98_map"
    @test sprint(Base.show, dustmap) == "SFD98Map($mapdir)"
    @test_throws ErrorException SFD98Map(randstring(10))
end

@testset "CSFD dust maps" begin
    # Reference values obtained from the dustmaps Python package (Green 2018)
    # using CSFDQuery at the following Galactic coordinates (l, b in radians):

    ref_l = [0.0, 1.0, 2.0, 3.0, 4.0, 0.0]
    ref_b = [0.5, 0.5, -0.5, 1.0, -1.0, π / 2]

    refebv = Float32[0.38901609, 0.04069592, 0.05450886, 0.01245212, 0.01313164, 0.01169674]

    # Reference mask values (bitmask from mask.fits)
    refmask = Int16[1, 5, 5, 5, 5, 1]

    dustmap = CSFDMap()

    for i in eachindex(refebv)
        @test dustmap(ref_l[i], ref_b[i]) ≈ refebv[i] rtol = 1e-5
    end
    # vectorized
    @test dustmap.(ref_l, ref_b) ≈ refebv rtol = 1e-5

    @testset "Mask dataproduct" begin
        maskmap = CSFDMap(dataproduct=:mask)
        for i in eachindex(refmask)
            @test maskmap(ref_l[i], ref_b[i]) == refmask[i]
        end
    end

    @testset "SkyCoords" begin
        for i in eachindex(refebv)
            @test dustmap(GalCoords(ref_l[i], ref_b[i])) ≈ refebv[i] rtol = 1e-5
            @test dustmap(convert(ICRSCoords, GalCoords(ref_l[i], ref_b[i]))) ≈ refebv[i] rtol = 1e-5
        end
    end

    @testset "Unitful" begin
        ref_l_u = ref_l * u"rad"
        ref_b_u = ref_b * u"rad"
        output = dustmap.(ref_l_u, ref_b_u)
        @test eltype(output) <: Gain
        @test ustrip.(output) ≈ refebv rtol = 1e-5
    end

    mapdir = datadep"csfd_map"
    @test sprint(Base.show, dustmap) == "CSFDMap($mapdir)"
    @test_throws ErrorException CSFDMap(randstring(10))
    @test_throws ErrorException CSFDMap(mapdir; dataproduct=:invalid)
end
