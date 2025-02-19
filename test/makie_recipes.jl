import Makie: Makie as M

law = CCM89()
aa = range(bounds(law)...; length=1_000)
xs = invum = map(DustExtinction.aa_to_invum, aa)
ys = map(law, aa)

@testset "PointBased" begin
    @test M.convert_arguments(M.PointBased(), law) == (M.Point{2, Float64}[[x, y] for (x, y) in zip(xs, ys)], )
    @test M.convert_arguments(M.PointBased(), aa,  law) == (M.Point{2, Float64}[[x, y] for (x, y) in zip(aa, ys)], )
    @test M.plottype(law) == M.Lines
    @test M.plottype(aa, law) == M.Lines
end

