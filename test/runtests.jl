using DustExtinction
using DustExtinction: bounds, checkbounds
using Test, Measurements, Unitful, UnitfulAstro, Random

Random.seed!(9994445781)

include("deprecate.jl")
include("color_laws.jl")
include("dust_maps.jl")
include("fittable_laws.jl")
include("mixture_laws.jl")

@testset "interfaces" begin
    for LAW in [CCM89, OD94, CAL00, GCC09, VCG04, FM90, G16, F99, F04, F19]
        @test bounds(LAW) == bounds(LAW())
        @test checkbounds(LAW, 1000) == checkbounds(LAW(), 1000)
        low, high = bounds(LAW)
        @test all(checkbounds.(LAW, rand(low:high, 1000)))
    end
    @test bounds(DustExtinction.ExtinctionLaw) == (0, Inf)
end

@testset "redden/deredden" begin
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

    # Av=0.3, Rv=3.1, law=ccm89
    ref_values = [0.6376288197244566, 0.6937993158943373, 0.7556003129054399, 0.7959216081971665, 0.8508142213648981, 0.6517368195297566, 0.6930246485833793, 0.7596467718783736, 0.7990835963390542, 0.8493553321733043, 0.9264835715958143, 0.9524307890106202, 0.9694546832913796, 0.7467375901587614, 0.8026674093639407, 0.847310247751628, 0.8771971841412943, 0.8125473475954554, 0.6442645437827309, 0.7127773090414787, 0.6749407318003154, 0.7249166081334651, 0.7578933014813821, 0.6461535564599767, 0.7288900774623704, 0.7900092679188825, 0.8409405763429691, 0.881972061390043, 0.6091986958411011, 0.7151404275193919, 0.7587309863781261, 0.7798780926274524, 0.8151697825084553, 0.8500202342098636, 0.7219200730427752, 0.7980286586611023, 0.8538476253948861]

    flux = ones(length(wave))
    output = @inferred broadcast((l, w, f)->redden(l, w, f; Av = 0.3), CCM89, wave, flux)
    @test output ≈ ref_values
    @test @inferred(broadcast((l, w, f)->deredden(l, w, f; Av = 0.3), CCM89, wave, output)) ≈ flux
    flux .= redden.(CCM89, wave, flux, Av = 0.3)
    @test flux ≈ ref_values


    # interfaces
    @test redden(CCM89, wave[1], flux[1]) == redden(CCM89(), wave[1], flux[1])
    @test redden(CCM89, wave[1], flux[1], Rv = 2.8) == redden(CCM89(Rv = 2.8), wave[1], flux[1])
    @test redden.(CCM89, wave, flux) == redden.(CCM89(), wave, flux)
    @test redden.(CCM89, wave, flux, Rv = 2.8) == redden.(CCM89(Rv = 2.8), wave, flux)

    # Measurements
    flux = ones(length(wave)) .± 0.1
    output = @inferred broadcast((l, w, f)->redden(l, w, f; Av = 0.3), CCM89, wave, flux)
    @test Measurements.value.(output) ≈ ref_values
    @test @inferred(broadcast((l, w, f)->deredden(l, w, f; Av = 0.3), CCM89, wave, output)) ≈ flux
    flux .= redden.(CCM89, wave, flux, Av = 0.3)
    @test flux ≈ ref_values

    # Unitful
    flux = ones(length(wave))u"Jy"
    wave = wave * u"angstrom"
    output = @inferred broadcast((l, w, f)->redden(l, w, f; Av = 0.3), CCM89, wave, flux)
    @test redden.(CCM89, wave, ustrip.(u"Jy", flux), Av = 0.3) ≈ ustrip.(u"Jy", output)
    @test ustrip.(output) ≈ ref_values
    @test @inferred(broadcast((l, w, f)->deredden(l, w, f; Av = 0.3), CCM89, wave, output)) ≈ flux
    @test deredden.(CCM89, wave, ustrip.(u"Jy", output), Av = 0.3) ≈ ustrip.(u"Jy", flux)
    flux .= redden.(CCM89, wave, flux, Av = 0.3)
    @test flux ≈ output
end
