import Makie: Makie as M

law = CCM89()
aa = range(bounds(law)...; length=1_000)
xs = invum = map(DustExtinction.aa_to_invum, aa)
ys = map(law, aa)
law_plot_args(wavs, vals) = (M.Point{2, Float64}[[x, y] for (x, y) in zip(wavs, vals)], )

@testset "PointBased" begin
    @test all(M.convert_arguments(M.PointBased(), law) .≈ law_plot_args(xs, ys))
    @test all(M.convert_arguments(M.PointBased(), aa, law) .≈ law_plot_args(aa, ys))
    @test M.plottype(law) == M.Lines
    @test M.plottype(aa, law) == M.Lines
end

lrange = range(-3, 3; length=400)
(lmin, lmax), llength, lstep = extrema(lrange), length(lrange) + 1, step(lrange)
lrange_plot = range(lmin - lstep/2, lmax + lstep/2; length=llength) |> collect

brange = range(-1, 1; length=300)
(bmin, bmax), blength, bstep = extrema(brange), length(brange) + 1, step(brange)
brange_plot = range(bmin - bstep/2, bmax + bstep/2; length=blength) |> collect

dustmap = SFD98Map()
m_plot = Float32[dustmap(li, bj) for li in lrange, bj in brange]
dustmap_plot_args = (lrange_plot, brange_plot, m_plot, )

@testset "CellGrid" begin
    @test all(M.convert_arguments(M.CellGrid(), dustmap) .≈ dustmap_plot_args)
    @test all(M.convert_arguments(M.CellGrid(), lrange,  brange, dustmap) .≈ dustmap_plot_args)
    @test M.plottype(dustmap) == M.Heatmap
    @test M.plottype(lrange, brange, dustmap) == M.Heatmap
end
