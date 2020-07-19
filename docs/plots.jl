using Plots, LaTeXStrings

import DustExtinction: ccm89_ca, ccm89_cb, od94_ca, od94_cb, cal00_invum, ccm89_invum, vcg04_invum, gcc09_invum, f99_invum, f04_invum, FM90, P92

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

#--------------------------------------------------------------------------------
# gcc09

w = range(3.3, 11.0, length=1000)
plot()
for rv in [2.0, 3.1, 4.0, 5.0, 6.0]
  m = gcc09_invum.(w, rv)
  plot!(w, m, label="Rv=$rv")
end
xlabel!(L"\mu m ^{-1}")
ylabel!("E(B-V)")
savefig(joinpath(dir, "gcc09_plot.svg"))

#--------------------------------------------------------------------------------
# vcg04

w = range(3.3, 8.0, length=1000)
plot()
for rv in [2.0, 3.1, 4.0, 5.0, 6.0]
  m = vcg04_invum.(w, rv)
  plot!(w, m, label="Rv=$rv")
end
xlabel!(L"\mu m ^{-1}")
ylabel!("E(B-V)")
savefig(joinpath(dir, "vcg04_plot.svg"))

#--------------------------------------------------------------------------------
# F04

w = range(0.3, 10.0, length=1000)
plot()
for rv in [2.0, 3.1, 4.0, 5.0, 6.0]
  m = f04_invum.(w, rv)
  plot!(w, m, label="Rv=$rv")
end
xlabel!(L"x\ \left[\mu m ^{-1}\right]")
ylabel!(L"A(x)/A(V)")
savefig(joinpath(dir, "F04_plot.svg"))

#--------------------------------------------------------------------------------
# FM90

w = range(3.8, 8.6, step = 0.001)
x = 1e4 ./ w
plot()

m1 = FM90().(x)
plot!(w, m1, label = "total")

m2 = FM90(c3=0.0, c4=0.0).(x)
plot!(w, m2, label = "linear term")

m3 = FM90(c1=0.0, c2=0.0, c4=0.0).(x)
plot!(w, m3, label = "bump term")

m4 = FM90(c1=0.0, c2=0.0, c3=0.0).(x)
plot!(w, m4, label = "FUV rise term")

xlabel!(L"\mu m ^{-1}")
ylabel!(L"E(\lambda - V)/E(B - V)")
savefig(joinpath(dir, "FM90_plot.svg"))

#--------------------------------------------------------------------------------
# P92

# plotting x in micrometers
# wave numbers generated in mu-1
w = exp10.(range(-3, 3, step = 0.001))

# wave generated in angstrom
x = 1e4 ./ w
# size was modified to accomodate the tally table
plot(size = (700, 500), xaxis = :log10, yaxis = :log10, xticks = ([0.001, 0.01, 0.1, 1, 10, 100, 1000]), yticks = ([0.001, 0.01, 0.1, 1, 10]))

m1 = P92().(x)
plot!(1 ./w, m1, label = "total")

m2 = P92(FUV_amp = 0, NUV_amp = 0, SIL1_amp = 0, SIL2_amp = 0, FIR_amp = 0).(x)
plot!(1 ./w, m2, label = "BKG only")

m3 = P92(NUV_amp=0.0, SIL1_amp=0.0, SIL2_amp=0.0, FIR_amp=0.0).(x)
plot!(1 ./w, m3, label = "BKG+FUV only")

m4 = P92(FUV_amp=0.0, SIL1_amp=0.0, SIL2_amp=0.0, FIR_amp=0.0).(x)
plot!(1 ./w, m4, label = "BKG+NUV only")

m5 = P92(FUV_amp = 0, NUV_amp = 0, SIL2_amp = 0).(x)
plot!(1 ./w, m5, label = "BKG+FIR+SIL1 only")

m6 = P92(FUV_amp=0., NUV_amp=0.0, SIL1_amp=0.0).(x)
plot!(1 ./w, m6, label = "BKG+FIR+SIL2 only")

m7 = P92(FUV_amp=0., NUV_amp=0.0, SIL1_amp=0.0, SIL2_amp=0.0).(x)
plot!(1 ./w, m7, label = "BKG+FIR only")

xlabel!(L"x\ \left[\mu m\right]")
ylabel!(L"A(X)/A(V)")
savefig(joinpath(dir, "P92_plot.svg"))

#--------------------------------------------------------------------------------
# F99

w = range(0.3, 10.0, length=1000)
plot()
for rv in [2.0, 3.1, 4.0, 5.0, 6.0]
  m = f99_invum.(w, rv)
  plot!(w, m, label="Rv=$rv")
end
xlabel!(L"x\ \left[\mu m ^{-1}\right]")
ylabel!(L"A(x)/A(V)")
savefig(joinpath(dir, "F99_plot.svg"))
