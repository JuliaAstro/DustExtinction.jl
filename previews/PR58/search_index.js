var documenterSearchIndex = {"docs":
[{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"using DustExtinction, CairoMakie\nusing DustExtinction: ExtinctionLaw\ninclude(\"./plotting.jl\")","category":"page"},{"location":"color_laws/#laws","page":"Color Laws","title":"Color laws","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"The following empirical laws allow us to model the reddening of light as it travels to us. The law you use should depend on the type of data you have and the goal of its use. CCM89 is very common for use in removing extinction from stellar observations, but CAL00, for instance, is suited for galaxies with massive stars. Look through the citations and documentation for each law to get a better idea of what sort of physics it targets.","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"DocTestSetup = quote\n    using DustExtinction, Random\n    Random.seed!(1)\n    ENV[\"UNITFUL_FANCY_EXPONENTS\"] = false\nend","category":"page"},{"location":"color_laws/#Usage","page":"Color Laws","title":"Usage","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"Color laws are constructed and then used as a function for passing wavelengths. Wavelengths are assumed to be in units of angstroms.","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"julia> CCM89(Rv=3.1)(4000)\n1.464555702942584","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"These laws can be applied across higher dimension arrays using the . operator","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"julia> CCM89(Rv=3.1).([4000, 5000])\n2-element Vector{Float64}:\n 1.464555702942584\n 1.1222468788993019","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"these laws return magnitudes, which we can apply directly to flux by mulitplication with a base-2.5 logarithmic system (because astronomers are fun):","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"f = f cdot 10 ^ -04A_vcdot mag","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"To make this easier, we provide a convenience redden and deredden functions for applying these color laws to flux measurements.","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"julia> wave = range(4000, 5000, length=4)\n4000.0:333.3333333333333:5000.0\n\njulia> flux = 1e-8 .* wave .+ 1e-2\n0.01004:3.3333333333333333e-6:0.01005\n\njulia> redden.(CCM89, wave, flux; Av=0.3)\n4-element Vector{Float64}:\n 0.00669864601545475\n 0.006918253926353551\n 0.007154659823737299\n 0.007370491272731541\n\njulia> deredden.(CCM89(Rv=3.1), wave, ans; Av=0.3) ≈ flux\ntrue\n","category":"page"},{"location":"color_laws/#Advanced-Usage","page":"Color Laws","title":"Advanced Usage","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"The color laws also have built-in support for uncertainties using Measurements.jl.","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"julia> using Measurements\n\njulia> CCM89(Rv=3.1).([4000. ± 10.5, 5000. ± 10.2])\n2-element Vector{Measurement{Float64}}:\n 1.4646 ± 0.0033\n 1.1222 ± 0.003\n","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"and also support units via Unitful.jl and its subsidiaries. Notice how the output type is now Unitful.Gain.","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"julia> using Unitful, UnitfulAstro\n\njulia> mags = CCM89(Rv=3.1).([4000u\"angstrom\", 0.5u\"μm\"])\n2-element Vector{Gain{Unitful.LogInfo{:Magnitude, 10, -2.5}, :?, Float64}}:\n 1.464555702942584 mag\n 1.1222468788993019 mag\n","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"You can even combine the two above to get some really nice workflows exploiting all Julia has to offer! This example shows how you could redden some OIR observational data with uncertainties in the flux density.","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"julia> using Measurements, Unitful, UnitfulAstro\n\njulia> wave = range(0.3, 1.0, length=5)u\"μm\"\n(0.3:0.175:1.0) μm\n\njulia> err = randn(length(wave))\n5-element Vector{Float64}:\n -0.07058313895389791\n  0.5314767537831963\n -0.806852326006714\n  2.456991333983293\n  1.1648740735275196\n\njulia> flux = @.(300 / ustrip(wave)^4 ± err)*u\"Jy\"\n5-element Vector{Quantity{Measurement{Float64}, 𝐌 𝐓^-2, Unitful.FreeUnits{(Jy,), 𝐌 𝐓^-2, nothing}}}:\n 37037.037 ± -0.071 Jy\n   5893.14 ± 0.53 Jy\n   1680.61 ± -0.81 Jy\n     647.6 ± 2.5 Jy\n     300.0 ± 1.2 Jy\n\njulia> redden.(CCM89, wave, flux; Av=0.3)\n5-element Vector{Quantity{Measurement{Float64}, 𝐌 𝐓^-2, Unitful.FreeUnits{(Jy,), 𝐌 𝐓^-2, nothing}}}:\n 22410.804 ± 0.043 Jy\n   4229.74 ± 0.38 Jy\n   1337.12 ± 0.64 Jy\n     554.3 ± 2.1 Jy\n     268.3 ± 1.0 Jy\n","category":"page"},{"location":"color_laws/#Parametric-Extinction-Laws","page":"Color Laws","title":"Parametric Extinction Laws","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"These laws are all parametrized by the selective extinction Rv. Mathematically, this is the ratio of the total extinction by the reddening","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"R_V = fracA_VE(B-V)","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"and is loosely associated with the size of the dust grains in the interstellar medium.","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"Index:","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"CCM89\nOD94\nCAL00\nVCG04\nGCC09\nF99\nF04\nF19","category":"page"},{"location":"color_laws/#Clayton,-Cardelli-and-Mathis-(1989)","page":"Color Laws","title":"Clayton, Cardelli and Mathis (1989)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"lplot(CCM89) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"CCM89","category":"page"},{"location":"color_laws/#DustExtinction.CCM89","page":"Color Laws","title":"DustExtinction.CCM89","text":"CCM89(;Rv=3.1)\n\nClayton, Cardelli and Mathis (1989) dust law.\n\nReturns E(B-V) in magnitudes at the given wavelength relative to the extinction at 5494.5 Å. The default support is [1000, 33333]. Outside of that range this will return 0. Rv is the selective extinction and is valid over [2, 6]. A typical value for the Milky Way is 3.1.\n\nReferences\n\nClayton,Cardelli and Mathis (1989)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#O'Donnell-1994","page":"Color Laws","title":"O'Donnell 1994","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"lplot(OD94) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"OD94","category":"page"},{"location":"color_laws/#DustExtinction.OD94","page":"Color Laws","title":"DustExtinction.OD94","text":"OD94(;Rv=3.1)\n\nO'Donnell (1994) dust law.\n\nThis is identical to the Clayton, Cardelli and Mathis (1989) dust law, except for different coefficients used in the optical (3030.3 Å to 9090.9 Å).\n\nReferences\n\nO'Donnell (1994)\n\nSee Also\n\nCCM89\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#Calzetti-et-al.-(2000)","page":"Color Laws","title":"Calzetti et al. (2000)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"lplot(CAL00) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"CAL00","category":"page"},{"location":"color_laws/#DustExtinction.CAL00","page":"Color Laws","title":"DustExtinction.CAL00","text":"CAL00(;Rv=4.05)\n\nCalzetti et al. (2000) Dust Law.\n\nReturns E(B-V) in magnitudes at the given wavelength. λ is the wavelength in Å and has support over [1200, 22000]. Outside of that range this will return 0.\n\nCalzetti et al. (2000) developed a recipe for dereddening the spectra of galaxies where massive stars dominate the radiation output. They found the best fit value for such galaxies was 4.05±0.80.\n\nReferences\n\nCalzetti et al. (2000)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#Valencic,-Clayton,-and-Gordon-(2004)","page":"Color Laws","title":"Valencic, Clayton, & Gordon (2004)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"lplot(VCG04) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"VCG04","category":"page"},{"location":"color_laws/#DustExtinction.VCG04","page":"Color Laws","title":"DustExtinction.VCG04","text":"VCG04(;Rv=3.1)\n\nValencic, Clayton, & Gordon (2004) dust law.\n\nThis model applies to the UV spectral region all the way to 912 Å. This model was not derived for the optical or NIR.\n\nReferences\n\nValencic, Clayton, & Gordon (2004)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#Gordon,-Cartledge,-and-Clayton-(2009)","page":"Color Laws","title":"Gordon, Cartledge, & Clayton (2009)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"lplot(GCC09) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"GCC09","category":"page"},{"location":"color_laws/#DustExtinction.GCC09","page":"Color Laws","title":"DustExtinction.GCC09","text":"GCC09(;Rv=3.1)\n\nGordon, Cartledge, & Clayton (2009) dust law.\n\nThis model applies to the UV spectral region all the way to 909.09 Å. This model was not derived for the optical or NIR.\n\nReferences\n\nGordon, Cartledge, & Clayton (2009)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#Fitzpatrick-(1999)","page":"Color Laws","title":"Fitzpatrick (1999)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"lplot(F99) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"F99","category":"page"},{"location":"color_laws/#DustExtinction.F99","page":"Color Laws","title":"DustExtinction.F99","text":"F99(;Rv=3.1)\n\nFitzpatrick (1999) dust law.\n\nReturns E(B-V) in magnitudes at the given wavelength relative to the extinction. This model applies to the UV and optical to NIR spectral range. The default support is [1000, 33333] Å. Outside of that range this will return 0. Rv is the selective extinction and is valid over [2, 6]. A typical value for the Milky Way is 3.1.\n\nReferences\n\nFitzpatrick (1999)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#Fitzpatrick-(2004)","page":"Color Laws","title":"Fitzpatrick (2004)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"lplot(F04) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"F04","category":"page"},{"location":"color_laws/#DustExtinction.F04","page":"Color Laws","title":"DustExtinction.F04","text":"F04(;Rv=3.1)\n\nFitzpatrick (2004) dust law.\n\nReturns E(B-V) in magnitudes at the given wavelength relative to the extinction. This model applies to the UV and optical to NIR spectral range. The default support is [1000, 33333] Å. Outside of that range this will return 0. Rv is the selective extinction and is valid over [2, 6]. A typical value for the Milky Way is 3.1.\n\nEquivalent to the F99 model with an updated NIR Rv dependence\n\nSee also Fitzpatrick & Massa (2007, ApJ, 663, 320)\n\nReferences\n\nFitzpatrick (2004)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#Fitzpatrick-(2019)","page":"Color Laws","title":"Fitzpatrick (2019)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"lplot(F19) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"F19","category":"page"},{"location":"color_laws/#DustExtinction.F19","page":"Color Laws","title":"DustExtinction.F19","text":"F19(;Rv=3.1)\n\nFitzpatrick (2019) dust law.\n\nReturns E(B-V) in magnitudes at the given wavelength relative to the extinction. This model applies to the UV and optical to NIR spectral range. The default support is [1149, 33333] Å. Outside of that range this will return 0. Rv is the selective extinction and is valid over [2, 6]. A typical value for the Milky Way is 3.1.\n\nFitzpatrick, Massa, Gordon et al. (2019, ApJ, 886, 108) model. Based on a sample of stars observed spectroscopically in the optical with HST/STIS.\n\nReferences\n\nFitzpatrick (2019)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#Maiz-Apellaniz-et-al.-(2014)","page":"Color Laws","title":"Maiz Apellaniz et al. (2014)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"lplot(M14) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"M14","category":"page"},{"location":"color_laws/#DustExtinction.M14","page":"Color Laws","title":"DustExtinction.M14","text":"M14(;Rv=3.1)\n\nMaiz Apellaniz et al (2014) Milky Way & LMC R(V) dependent model.\n\nReturns E(B-V) in magnitudes at the given wavelength relative to the extinction. The published UV extinction curve is identical to Clayton, Cardelli, and Mathis (1989, CCM). Forcing the optical section to match smoothly with CCM introduces a non-physical feature at high values of R5495 around 3.9 inverse microns; see section 5 in Maiz Apellaniz et al. (2014) for more discussion. For that reason, we provide the M14 model only through 3.3 inverse microns, the limit of the optical in CCM. Outside of that range this will return 0. Rv is the selective extinction and is valid over [2, 6]. A typical value for the Milky Way is 3.1. R5495 = A(5485)/E(4405-5495) Spectral equivalent to photometric R(V).\n\nReferences\n\nMaiz Apellaniz et al. (2014)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#API/Reference","page":"Color Laws","title":"API/Reference","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"redden\nderedden\nDustExtinction.ExtinctionLaw\nDustExtinction.bounds\nDustExtinction.checkbounds","category":"page"},{"location":"color_laws/#DustExtinction.redden","page":"Color Laws","title":"DustExtinction.redden","text":"redden(::ExtinctionLaw, wave, flux; Av=1)\nredden(::Type{ExtinctionLaw}, wave, flux; Av=1, law_kwargs...)\n\nRedden the given flux using the given extinction law at the given wavelength.\n\nIf wave is <:Real then it is expected to be in angstrom and if it is <:Unitful.Quantity it will be automatically converted. Av is the total extinction value. The extinction law can be a constructed struct or a Type. If it is a Type, law_kwargs will be passed to the constructor.\n\nExamples\n\njulia> wave = 3000; flux = 1000;\n\njulia> redden(CCM89, wave, flux; Rv=3.1)\n187.38607779757183\n\njulia> redden(CCM89(Rv=3.1), wave, flux; Av=2)\n35.11354215235764\n\nSee Also\n\nderedden\n\n\n\n\n\n","category":"function"},{"location":"color_laws/#DustExtinction.deredden","page":"Color Laws","title":"DustExtinction.deredden","text":"deredden(::ExtinctionLaw, wave, flux; Av=1)\nderedden(::Type{ExtinctionLaw}, wave, flux; Av=1, law_kwargs...)\n\nDeredden the given flux using the given extinction law at the given wavelength.\n\nIf wave is <:Real then it is expected to be in angstrom and if it is <:Unitful.Quantity it will be automatically converted. Av is the total extinction value. The extinction law can be a constructed struct or a Type. If it is a Type, law_kwargs will be passed to the constructor.\n\nExamples\n\njulia> wave = 3000; flux = 187.386;\n\njulia> deredden(CCM89, wave, flux; Rv=3.1)\n999.9995848273642\n\njulia> deredden(CCM89(Rv=3.1), wave, flux; Av=2)\n5336.573541539394\n\nSee Also\n\nredden\n\n\n\n\n\n","category":"function"},{"location":"color_laws/#DustExtinction.ExtinctionLaw","page":"Color Laws","title":"DustExtinction.ExtinctionLaw","text":"abstract type DustExtinction.ExtinctionLaw\n\nThe abstract supertype for dust extinction laws. See the extended help (??DustExtinction.ExtinctionLaw from the REPL) for more information about the interface.\n\nExtended Help\n\nInterface\n\nHere's how to make a new extinction law, called MyLaw\n\nCreate your struct. We strongly recommend using Parameters.jl to facilitate creating keyword argument constructors if your model is parameterized, which allows convenient usage with redden and deredden.\nstruct MyLaw <: DustExtinction.ExtinctionLaw end\n(Optional) Define the limits. This will default to (0, Inf). Currently, this is used within the DustExtinction.checkbounds function and in the future will be used for plotting recipes.\nDustExtinction.bounds(::Type{<:MyLaw}) = (min, max)\nDefine the law. You only need to provide one function which takes wavelength as angstrom. If your law is naturally written for inverse-micron, there is a helper function aa_to_invum.\n(::MyLaw)(wavelength::Real)\n(Optional) enable Unitful.jl support by adding this function. If you are building a new law within DustExtinction.jl you can add your law to the code-gen list inside DustExtinction.jl/src/DustExtinction.jl.\n(l::MyLaw)(wavelength::Unitful.Quantity) = l(ustrip(u\"angstrom\", wavelength)) * u\"mag\"\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#DustExtinction.bounds","page":"Color Laws","title":"DustExtinction.bounds","text":"DustExtinction.bounds(::ExtinctionLaw)::Tuple\nDustExtinction.bounds(::Type{<:ExtinctionLaw})::Tuple\n\nGet the natural wavelengths bounds for the extinction law, in angstrom\n\n\n\n\n\n","category":"function"},{"location":"color_laws/#DustExtinction.checkbounds","page":"Color Laws","title":"DustExtinction.checkbounds","text":"DustExtinction.checkbounds(::ExtinctionLaw, wavelength)::Bool\nDustExtinction.checkbounds(::Type{<:ExtinctionLaw, wavelength}::Bool\n\nHelper function that uses DustExtinction.bounds to return whether the given wavelength is in the support for the law.\n\n\n\n\n\n","category":"function"},{"location":"color_laws/#Fittable-Extinction-Laws","page":"Color Laws","title":"Fittable Extinction Laws","text":"","category":"section"},{"location":"color_laws/#Fitzpatrick-and-Massa-(1990)","page":"Color Laws","title":"Fitzpatrick & Massa (1990)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"lplot(FM90) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"FM90","category":"page"},{"location":"color_laws/#DustExtinction.FM90","page":"Color Laws","title":"DustExtinction.FM90","text":"FM90(;c1=0.10, c2=0.70, c3=3.23, c4=0.41, x0=4.60, gamma=0.9)\nFM90(coeffs, x0=4.60, gamma=0.9)\n\nFitzpatrick & Massa (1990) generative model for ultraviolet dust extinction. The default values are the published values for the Milky Way average.\n\nParameters\n\nc1 - y-intercept of linear term\nc2 - slope of liner term\nc3 - amplitude of 2175 Å bump\nc4 - amplitude of FUV rise\nx0 - centroid of 2175 Å bump\ngamma - width of 2175 Å bump\n\nIf λ is a Unitful.Quantity it will be automatically converted to Å and the returned value will be UnitfulAstro.mag.\n\nExamples\n\njulia> model = FM90(c1=0.2, c2=0.7, c3=3.23, c4=0.41, x0=4.6, gamma=0.99);\n\njulia> model(1500)\n5.2521258452800135\n\njulia> FM90()(1500)\n5.152125845280013\n\njulia> FM90(c1=0.2, c2=0.7, c3=3.23, c4=0.41, x0=4.6, gamma=0.99).([1000, 1200, 1800])\n3-element Vector{Float64}:\n 12.562237969522851\n  7.769215017329513\n  4.890128210972148\n\n\nExtended Help\n\nThe model has form c_1 + c_2x + c_3D(x γ x_0) + c_4 F(x) where x is the wavenumber in inverse microns, D(x) is a Drude profile (modified Lorentzian) used to model the 2175 Å bump with the scale-free parameters x_0 (central wavenumber) and γ (damping coefficient), and F(x), a piecewise function for the far-UV. Note that the coefficients will change the overall normalization, possibly changing the expected behavior of reddening via the parameter A_V.\n\nReferences\n\nFitzpatrick & Massa (1990)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#Mixture-Extinction-Laws","page":"Color Laws","title":"Mixture Extinction Laws","text":"","category":"section"},{"location":"color_laws/#Gordon-et-al.-(2003)","page":"Color Laws","title":"Gordon et al. (2003)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"G03_SMCBar\nG03_LMCAve","category":"page"},{"location":"color_laws/#DustExtinction.G03_SMCBar","page":"Color Laws","title":"DustExtinction.G03_SMCBar","text":"G03_SMCBar(;Rv=2.74) <Internal function>\n\nGordon et al. (2003) SMCBar Average Extinction Curve.\n\nThe observed A(lambda)/A(V) values at 2.198 and 1.25 microns were changed to provide smooth interpolation as noted in Gordon et al. (2016, ApJ, 826, 104)\n\nReference\n\nGordon et al. (2003)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#DustExtinction.G03_LMCAve","page":"Color Laws","title":"DustExtinction.G03_LMCAve","text":"G03_LMCAve(;Rv=3.41) <Internal function>\n\nGordon et al. (2003) LMCAve Average Extinction Curve.\n\nReference\n\nGordon et al. (2003)\n\n\n\n\n\n","category":"type"},{"location":"color_laws/#Gordon-et-al.-(2016)","page":"Color Laws","title":"Gordon et al. (2016)","text":"","category":"section"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"mplot(G16, (2.0, 3.1, 4.0, 5.0, 6.0), 1.0) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"mplot(G16, 3.1, 0.0:0.2:1.0) # hide","category":"page"},{"location":"color_laws/","page":"Color Laws","title":"Color Laws","text":"G16","category":"page"},{"location":"color_laws/#DustExtinction.G16","page":"Color Laws","title":"DustExtinction.G16","text":"G16(;Rv=3.1, f_A=1.0)\n\nGordon et al. (2016) Milky Way, LMC, & SMC R(V) and f_A dependent model\n\nReturns E(B-V) in magnitudes at the given wavelength relative to the extinction. This is mixture model between the F99 R(V) dependent model (component A) and the G03_SMCBar model (component B). The default support is [1000, 33333] Å. Outside of that range this will return 0. Rv is the selective extinction and is valid over [2, 6]. A typical value for the Milky Way is 3.1.\n\nReferences\n\nGordon et al. (2016)\n\n\n\n\n\n","category":"type"},{"location":"plotting/#plotting","page":"Plotting","title":"Plotting","text":"","category":"section"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"DustExtinction.jl provides basic functionality for creating plots with Makie.jl, which can then be composed together to form more complex figures. Below are a few examples using the following primitives from Makie.jl: lines, scatter, and heatmap, along with their mutating equivalents.","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"note: Note\nBy default, all plots adopt the axis limits defined by DustExtinction.bounds, but this is easily modified in the Makie figures.","category":"page"},{"location":"plotting/#Model-plot-example","page":"Plotting","title":"Model plot example","text":"","category":"section"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"Makie's line and scatter plots work out-of-the box for all color laws, which can then be used to build up more complex plots.","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"using DustExtinction, CairoMakie\nusing DustExtinction: aa_to_invum\nusing Measurements","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"wavs = let\n    x = range(2_000, 3_000; length=1_000)\n    x .± 1e5 * inv.(x)\nend # Å","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"invum = aa_to_invum.(wavs) # μm","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"model = CCM89()","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"lines(model)","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"lines(model; axis=(; limits=(3, 6, 0, 5)))","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"extinction = model.(wavs) # mag","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"lines(wavs, extinction)","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"wavs_sampled, extinction_sampled = let\n    N_samples = 7\n    wavs[range(begin, step=end ÷ N_samples; length=N_samples)],\n    extinction[range(begin, step=end ÷ N_samples; length=N_samples)]\nend","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"fig, ax, p = band(wavs, extinction; alpha=0.5, label=\"model uncertainty\")\n\nlines!(ax, wavs, extinction; label=\"model: CCM89\")\n\nscatter!(ax, wavs_sampled, extinction_sampled; color=:orange, label=\"observations\")\n\n# Currently ambiguous for both x and y being Measurements\nerrorbars!(ax, Measurements.value.(wavs_sampled), extinction_sampled;\n    whiskerwidth = 10,\n    color = :orange,\n    label = \"obs. uncertainty\",\n)\n\naxislegend()\nax.xlabel = \"Wavelength [Å]\"\nax.ylabel = \"A(x) / A(V) [mag]\"\n\nfig","category":"page"},{"location":"plotting/#Dust-map-example","page":"Plotting","title":"Dust map example","text":"","category":"section"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"# Example here\n\n\n\n\n2 + 2","category":"page"},{"location":"plotting/","page":"Plotting","title":"Plotting","text":"tip: Tip\nSee plotting.jl for more plotting examples. The convenience functions defined there are used to generate the other figures shown in this documentation.","category":"page"},{"location":"#DustExtinction.jl","page":"Home","title":"DustExtinction.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: GitHub) (Image: Build Status) (Image: PkgEval) (Image: Coverage) (Image: License)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package implements different empirical dust measurements for use in astronomy. This package is written in pure Julia and is built with first-class support with Unitful.jl and Measurements.jl.","category":"page"},{"location":"#About","page":"Home","title":"About","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Extinction describes the effect of dust grains on observations of stars in space. Light that travels through dust is absorbed and scatterred as natural processes of light's interactions with materials. This obfuscation can be modeled and removed from our data in order to more properly retrieve the star's flux. When dealing with multiple stars, or large clusters or galaxies, this process is considered dust attenuation and is not provided for in this package.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"From the REPL, press ] to enter Pkg mode","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add DustExtinction\n\njulia> using DustExtinction","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using DustExtinction\n\njulia> CCM89(Rv=3.1)(4000)\n1.464555702942584","category":"page"},{"location":"","page":"Home","title":"Home","text":"For more examples, view the Color Laws and Dust Maps sections.","category":"page"},{"location":"#Citations","page":"Home","title":"Citations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are various citations relevant to this work. Please be considerate when using this work or any derivate of it by adding the appropriate citations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Law Reference BibTeX\nCCM89 Clayton, Cardelli and Mathis (1989) download\nOD94 O'Donnell (1994) download\nCAL00 Calzetti et al. (2000) download\nVCG04 Valencic, Clayton, & Gordon (2004) download\nGCC09 Gordon, Cartledge, & Clayton (2009) download\nFM90 Fitzpatrick & Massa (1990) download\nG16 Gordon et al (2016)  download\nSFD98Map Schlegel, Finkbeiner and Davis (1998) download\nF99 Fitzpatrick (1999) download\nF04 Fitzpatrick (2004) download\nF19 Fitzpatrick (2019) download\nM14 Maiz Apellaniz et al. (2014) download","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [DustExtinction]\nOrder = [:function, :type]","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you are interested in contributing, feel free to make a pull request or open an issue for discussion.","category":"page"},{"location":"dust_maps/","page":"Dust Maps","title":"Dust Maps","text":"using DustExtinction, CairoMakie\nusing DustExtinction: ExtinctionLaw\ninclude(\"./plotting.jl\")","category":"page"},{"location":"dust_maps/#maps","page":"Dust Maps","title":"Dust Maps","text":"","category":"section"},{"location":"dust_maps/#Usage","page":"Dust Maps","title":"Usage","text":"","category":"section"},{"location":"dust_maps/","page":"Dust Maps","title":"Dust Maps","text":"DocTestSetup = quote\n    using DustExtinction, Random\n    Random.seed!(1)\nend","category":"page"},{"location":"dust_maps/","page":"Dust Maps","title":"Dust Maps","text":"julia> dustmap = SFD98Map();\n\njulia> dustmap(0, 2)\n0.020303287464050277\n\njulia> l = range(-π, π, length=5)\n-3.141592653589793:1.5707963267948966:3.141592653589793\n\njulia> b = range(-π/2, π/2, length=5)\n-1.5707963267948966:0.7853981633974483:1.5707963267948966\n\njulia> [dustmap(l[i], b[j]) for i in 1:length(l), j in 1:length(b)]\n5×5 Matrix{Float64}:\n 0.0159853  0.105782    1.40486  0.0158918  0.0119615\n 0.0159853  0.0268289   3.47788  0.0654852  0.0119615\n 0.0159853  0.0343457  99.6976   0.103875   0.0119615\n 0.0159853  0.0432165   2.60569  0.0178195  0.0119615\n 0.0159853  0.105782    1.40486  0.0158918  0.0119615\n","category":"page"},{"location":"dust_maps/","page":"Dust Maps","title":"Dust Maps","text":"dplot() # hide","category":"page"},{"location":"dust_maps/#Advanced-Usage","page":"Dust Maps","title":"Advanced Usage","text":"","category":"section"},{"location":"dust_maps/","page":"Dust Maps","title":"Dust Maps","text":"Our dust maps also have native support for Unitful.jl and Measurements.jl.","category":"page"},{"location":"dust_maps/","page":"Dust Maps","title":"Dust Maps","text":"julia> using Measurements, Unitful\n\njulia> l = 45u\"°\"; b=0u\"°\";\n\njulia> dustmap = SFD98Map();\n\njulia> dustmap(l, b)\n6.4290331211742355 mag\n\njulia> l = l ± 0.1u\"°\"; b = b ± 0.3u\"°\";\n\njulia> dustmap(l, b)\n6.4 ± 5.7 mag\n","category":"page"},{"location":"dust_maps/#API/Reference","page":"Dust Maps","title":"API/Reference","text":"","category":"section"},{"location":"dust_maps/","page":"Dust Maps","title":"Dust Maps","text":"SFD98Map","category":"page"},{"location":"dust_maps/#DustExtinction.SFD98Map","page":"Dust Maps","title":"DustExtinction.SFD98Map","text":"SFD98Map([mapdir])\n\nSchlegel, Finkbeiner and Davis (1998) dust map.\n\nThe first time this is constructed, the data files required will be downloaded and stored in a directory following the semantics of DataDeps.jl. To avoid being asked to download the files, set the environment variable DATADEPS_ALWAYS_ACCEPT to true. You can also provide the directory of the two requisite files manually instead of relying on DataDeps.jl. Internally, this type keeps the FITS files defining the map open, speeding up repeated queries for E(B-V) values.\n\nReferences\n\nSchlegel, Finkbeiner and Davis (1998)\n\n\n\n\n\n(dustmap::SFD98Map)(l::Real, b::Real)\n(dustmap::SFD98Map)(l::Quantity, b::Quantity)\n\nGet E(B-V) value from a SFD98Map instance at galactic coordinates (l, b), given in radians. Uses bilinear interpolation between pixel values. If l and b are Unitful.Quantity they will be converted to radians and the output will be given as UnitfulAstro.mag.\n\nExample\n\njulia> using DustExtinction\n\njulia> m = SFD98Map();\n\njulia> m(1, 2)\n0.013439524544325624\n\njulia> l = 0:0.5:2; b = 0:0.5:2;\n\njulia> m.(l, b)\n5-element Vector{Float64}:\n 99.69757461547852\n  0.10180447359074371\n  0.019595484241066132\n  0.010238757633890877\n  0.01862100327420125\n\n\n\n\n\n","category":"type"}]
}
