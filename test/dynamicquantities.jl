import DynamicQuantities as DQ

# NOTE: `DynamicQuantities` and `Unitful` both export `@u_str`/`ustrip`, and both
# are loaded into this test module, so we fully qualify the DynamicQuantities
# versions (e.g. `DQ.u"angstrom"`) to avoid the name collision.
@testset "DynamicQuantities" begin
    law = CCM89(Rv = 3.1)

    @testset "law evaluation" begin
        wave = 3000.0
        # A wavelength quantity is converted to angstrom; result is plain
        # magnitudes (DynamicQuantities has no logarithmic `mag` unit).
        @test law(wave * DQ.u"angstrom") ≈ law(wave)
        @test law(300 * DQ.u"nm") ≈ law(wave)  # 3000 Å == 300 nm
        @test law(wave * DQ.u"angstrom") isa Real
    end

    @testset "redden/deredden" begin
        wave = 3000.0
        flux = 1000.0
        wave_q = wave * DQ.u"angstrom"

        # Real flux, unitful wavelength
        @test redden(law, wave_q, flux; Av = 0.3) ≈ redden(law, wave, flux; Av = 0.3)
        @test deredden(law, wave_q, flux; Av = 0.3) ≈ deredden(law, wave, flux; Av = 0.3)
        @test redden(CCM89, wave_q, flux; Rv = 3.1) ≈ redden(CCM89, wave, flux; Rv = 3.1)

        # Unitful flux preserves its units and round-trips
        flux_q = flux * DQ.u"W"
        reddened = redden(law, wave_q, flux_q; Av = 0.3)
        @test DQ.dimension(reddened) == DQ.dimension(flux_q)
        @test DQ.ustrip(reddened) ≈ redden(law, wave, flux; Av = 0.3)
        @test deredden(law, wave_q, reddened; Av = 0.3) ≈ flux_q

        # in-place with unitful flux round-trips
        flux_mut = [flux_q, flux_q]
        redden!(law, [wave_q, wave_q], flux_mut; Av = 0.3)
        @test flux_mut ≈ [reddened, reddened]
        deredden!(law, [wave_q, wave_q], flux_mut; Av = 0.3)
        @test flux_mut ≈ [flux_q, flux_q]
    end
end
