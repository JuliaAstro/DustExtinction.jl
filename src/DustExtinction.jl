module DustExtinction

export ccm89,
       od94

# Optical coefficients
const ccm89_ca = [1. 0.17699 -0.50447 -0.02427 0.72085 0.01979 -0.77530 0.32999]
const ccm89_cb = [0. 1.41338 2.28305 1.07233 -5.38434 -0.62251 5.30260 -2.09002]
const od94_ca = [1. 0.104 -0.609  0.701  1.137 -1.718 -0.827   1.647 -0.505]
const od94_cb = [0. 1.952  2.908 -3.989 -7.985 11.102  5.491 -10.805  3.347]

function ccm89like(w::FloatingPoint, r_v, c_a, c_b)
    x = 1.e4 / w
    a = 0.
    b = 0.
    if x < 0.3
        error("wavelength out of range")
    elseif x < 1.  # Near IR
        y = x^1.61
        a = 0.574 * y
        b = -0.527 * y
    elseif x < 3.3  # Optical
        y = x - 1.82
        yn = 1.
        a = c_a[1]
        b = c_b[1]
        for i = 2:length(c_a)
            yn *= y
            a += c_a[i] * yn
            b += c_b[i] * yn
        end
    elseif x < 8.  # UV
        a =  1.752 - 0.316*x - (0.104 / ((x-4.67)^2 + 0.341))
        b = -3.090 + 1.825*x + (1.206 / ((x-4.62)^2 + 0.263))
        if x > 5.9
            y = x - 5.9
            y2 = y * y
            y3 = y2 * y
            a += -0.04473*y2 - 0.009779*y3
            b += 0.2130*y2 + 0.1207*y3
        end
    elseif x < 11.
        y = x - 8.
        y2 = y * y
        y3 = y2 * y
        a = -0.070*y3 + 0.137*y2 - 0.628*y - 1.073
        b = 0.374*y3 - 0.420*y2 + 4.257*y + 13.670
    else
        error("wavelength out of range")
    end

    a + b/r_v
end

ccm89(w::FloatingPoint, r_v=3.1) = ccm89like(w, r_v, ccm89_ca, ccm89_cb)
od94(w::FloatingPoint, r_v=3.1) = ccm89like(w, r_v, od94_ca, od94_cb)

for f = (:ccm89, :od94)
    @eval begin
        function ($f){T<:FloatingPoint}(w::AbstractArray{T,1}, r_v=3.1)
            [ ($f)(w[i], r_v) for i=1:length(w) ]
        end
    end
end

end # module
