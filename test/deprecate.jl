
@testset "v0.4 deprecations" begin
    x_inv_microns = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.6, 4.0, 3.0, 2.78, 2.27, 1.82, 1.43, 1.11, 0.8, 0.63, 0.46]
    wave = 1e4 ./ x_inv_microns
    for law in [ccm89, od94]
        @test_deprecated law(wave)
    end
end

@testset "v0.7 deprecations" begin
    wave = 3000
    flux = 1

    @test_deprecated redden(wave, flux, 0.3)
    @test_deprecated deredden(wave, flux, 0.3)

    for law in [ccm89, od94, cal00, vcg04, gcc09]
        @test_deprecated law(wave)
    end
end
