using Plots, LaTeXStrings
import DustExtinction: ccm89_ca, ccm89_cb, od94_ca, od94_cb, cal00_invum, ccm89_invum

dir = joinpath(@__DIR__, "src", "assets")

#--------------------------------------------------------------------------------
# ccm89 

w = range(0.3, 10.0, length=1000)
plot()
for rv in [2.0, 3.1, 4.0, 5.0, 6.0]
  m = ccm89_invum.(w, rv, Ref(ccm89_ca), Ref(ccm89_cb))
  plot!(w, m, label="Rv=$rv")
end
xlabel!(L"\mu m ^{-1}")
ylabel!("E(B-V)")
savefig(joinpath(dir, "ccm89_plot.svg"))

#--------------------------------------------------------------------------------
# od94 

w = range(0.3, 10.0, length=1000)
plot()
for rv in [2.0, 3.1, 4.0, 5.0, 6.0]
  m = ccm89_invum.(w, rv, Ref(od94_ca), Ref(od94_cb))
  plot!(w, m, label="Rv=$rv")
end
xlabel!(L"\mu m ^{-1}")
ylabel!("E(B-V)")
savefig(joinpath(dir, "od94_plot.svg"))

#--------------------------------------------------------------------------------
# cal00 

w = range(0.46, 8.3, length=1000)
plot()
for rv in [2.0, 3.0, 4.05, 5.0, 6.0]
  m = cal00_invum.(w, rv)
  plot!(w, m, label="Rv=$rv")
end
xlabel!(L"\mu m ^{-1}")
ylabel!("E(B-V)")
savefig(joinpath(dir, "cal00_plot.svg"))

#--------------------------------------------------------------------------------
# sfd98 

l = range(-pi, pi, length=400)
b = range(-pi/64, pi/64, length=300)
dustmap = SFD98Map()
m = [dustmap(li, bj) for li in l, bj in b]
heatmap(l, b, m, label="", transpose=true, colorbar_title="E(B-V)")
xlabel!("l (rad)")
ylabel!("b (rad)")
savefig(joinpath(dir, "sfd98_plot.svg"))
