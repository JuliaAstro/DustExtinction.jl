# DustExtinction.jl

This package implements different empirical dust measurements for use in astronomy. This package is written in pure Julia and is built with first-class support with [`Unitful.jl`](https://github.com/painterqubits/unitful.jl) and [`Measurements.jl`](https://github.com/juliaphysics/measurements.jl).

## Installation

From the REPL, press `]` to enter Pkg mode

```julia
(v 1.1) pkg> add DustExtinction
[...]

julia> using DustExtinction
```

## Usage

```jldoctest
julia> using DustExtinction

julia> ccm89(4000., 3.1)
1.4645557029425842

```

For more examples, view the [Color Laws](@ref laws) and [Dust Maps](@ref maps) sections. 

## API/Reference

```@index
Modules = [DustExtinction]
Order = [:function, :type]
```

## Citations

There are various citations relevant to this work. Please be considerate when using this work or any derivate of it by adding the appropriate citations.

[`ccm89`](@ref) - 
[Cardelli, Clayton and Mathis (1989)](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C/abstract)

```
@ARTICLE{1989ApJ...345..245C,
       author = {{Cardelli}, Jason A. and {Clayton}, Geoffrey C. and {Mathis}, John S.},
        title = "{The Relationship between Infrared, Optical, and Ultraviolet Extinction}",
      journal = {\apj},
     keywords = {Infrared Spectra, Interstellar Extinction, Ultraviolet Spectra, Visible Spectrum, Computational Astrophysics, Interstellar Matter, Iue, Astrophysics, INTERSTELLAR: MATTER, ULTRAVIOLET: SPECTRA},
         year = "1989",
        month = "Oct",
       volume = {345},
        pages = {245},
          doi = {10.1086/167900},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

[`od94`](@ref) - [O'Donnell (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...422..158O/abstract)
```
@ARTICLE{1994ApJ...422..158O,
       author = {{O'Donnell}, James E.},
        title = "{R v-dependent Optical and Near-Ultraviolet Extinction}",
      journal = {\apj},
     keywords = {Interstellar Extinction, Light (Visible Radiation), Near Infrared Radiation, Ultraviolet Radiation, Astronomical Photometry, Iue, Astrophysics, ISM: DUST, EXTINCTION},
         year = "1994",
        month = "Feb",
       volume = {422},
        pages = {158},
          doi = {10.1086/173713},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1994ApJ...422..158O},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

[`cal00`](@ref) - [Calzetti et al. (2000)](https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C/abstract)
```
@ARTICLE{2000ApJ...533..682C,
       author = {{Calzetti}, Daniela and {Armus}, Lee and {Bohlin}, Ralph C. and
         {Kinney}, Anne L. and {Koornneef}, Jan and {Storchi-Bergmann}, Thaisa},
        title = "{The Dust Content and Opacity of Actively Star-forming Galaxies}",
      journal = {\apj},
     keywords = {GALAXIES: STARBURST, INFRARED: GALAXIES, INFRARED: ISM: CONTINUUM, ISM: DUST, EXTINCTION, Astrophysics},
         year = "2000",
        month = "Apr",
       volume = {533},
       number = {2},
        pages = {682-695},
          doi = {10.1086/308692},
archivePrefix = {arXiv},
       eprint = {astro-ph/9911459},
 primaryClass = {astro-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

[`SFD98Map`](@ref) - [Schlegel, Finkbeiner and Davis (1998)](https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract)
```
@ARTICLE{1998ApJ...500..525S,
       author = {{Schlegel}, David J. and {Finkbeiner}, Douglas P. and {Davis}, Marc},
        title = "{Maps of Dust Infrared Emission for Use in Estimation of Reddening and Cosmic Microwave Background Radiation Foregrounds}",
      journal = {\apj},
     keywords = {COSMOLOGY: DIFFUSE RADIATION, COSMOLOGY: COSMIC MICROWAVE BACKGROUND, ISM: DUST, EXTINCTION, INTERPLANETARY MEDIUM, INFRARED: ISM: CONTINUUM, Cosmology: Cosmic Microwave Background, Cosmology: Diffuse Radiation, ISM: Dust, Extinction, Infrared: ISM: Continuum, Interplanetary Medium, Astrophysics},
         year = "1998",
        month = "Jun",
       volume = {500},
       number = {2},
        pages = {525-553},
          doi = {10.1086/305772},
archivePrefix = {arXiv},
       eprint = {astro-ph/9710327},
 primaryClass = {astro-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## Contributing

If you are interested in contributing, feel free to make a pull request or open an issue for discussion. 
