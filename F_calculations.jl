include("cosmology_vars.jl")
include("nu_osc_params.jl")

using DelimitedFiles
using Plots
using Interpolations
using LaTeXStrings
using QuadGK
using SpecialFunctions
using BenchmarkTools
# using Polynomials

function F0_tint_func(F0vec::Vector{Float64})
    es = range(0, 100, 2000)
    F0int_ne = Interpolations.interpolate((vec(es),), F0vec, Gridded(Linear()))
    return extrapolate(F0int_ne, 0.0)
end

F0s_vec = readdlm("F0s_vec.txt", comments=true)

F0_νe_270sm, F0_νebar_270sm, F0_νx_270sm = F0_tint_func(vec(F0s_vec[:, 1])), F0_tint_func(vec(F0s_vec[:, 2])), F0_tint_func(vec(F0s_vec[:, 3]))
F0_νe_112sm, F0_νebar_112sm, F0_νx_112sm = F0_tint_func(vec(F0s_vec[:, 4])), F0_tint_func(vec(F0s_vec[:, 5])), F0_tint_func(vec(F0s_vec[:, 6]))
F0_νe_bh, F0_νebar_bh, F0_νx_bh = F0_tint_func(vec(F0s_vec[:, 7])), F0_tint_func(vec(F0s_vec[:, 8])), F0_tint_func(vec(F0s_vec[:, 9]))

function F0(E, β, sm)
    if β == "e" && sm == "small"
        return F0_νe_112sm(E)
    elseif β == "e" && sm == "large"
        return F0_νe_270sm(E)
    elseif β == "e" && sm == "bh"
        return F0_νe_bh(E)
    elseif β == "ebar" && sm == "small"
        return F0_νebar_112sm(E)
    elseif β == "ebar" && sm == "large"
        return F0_νebar_270sm(E)
    elseif β == "ebar" && sm == "bh"
        return F0_νebar_bh(E)   
    elseif β == "x" && sm == "small"
        return F0_νx_112sm(E)
    elseif β == "x" && sm == "large"
        return F0_νx_270sm(E)
    elseif β == "x" && sm == "bh"
        return F0_νx_bh(E)
    else
        return 0
    end
end

# Oscillations thru the SN medium
# Accounting for oscillations thru the SN medium

s12 = 0.297
c12 = 1 - s12
PH = 0

# ordering = "NO" (normal ordering) or "IO" (inverted ordering)
function F(E, β, sm, ordering)
    if ordering == "NO"
        if β == "e"
            return F0(E, "x", sm)
        elseif β == "ebar"
            return c12*F0(E, "ebar", sm) + s12*F0(E, "x", sm)
        elseif β == "x"
            return 0.5*(F0(E, "e", sm) + F0(E, "x", sm))
        elseif β == "xbar"
            return 0.5*(s12*F0(E, "ebar", sm) + (1 + c12)*F0(E, "x", sm))
        else
            return 0
        end
    elseif ordering == "IO"
        if β == "e"
            return s12*F0(E, "e", sm) + c12*F0(E, "x", sm)
        elseif β == "ebar"
            return F0(E, "x", sm)
        elseif β == "x"
            return 0.5*(c12*F0(E, "e", sm) + (1 + s12)*F0(E, "x", sm))
        elseif β == "xbar"
            return 0.5*(F0(E, "ebar", sm) + F0(E, "x", sm))
        else
            return 0
        end
    else
        return 0
    end
end

# In the mass basis now: i = 1, 2, 3, nubar = true or false
function Fmass_old(E, i, sm, ordering, nubar)
    if nubar==false
        return Usqred(ordering)[1, i]*F(E, "e", sm, ordering) + (1 - Usqred(ordering)[1, i])*F(E, "x", sm, ordering)
    elseif nubar==true
        return Usqred(ordering)[1, i]*F(E, "ebar", sm, ordering) + (1 - Usqred(ordering)[1, i])*F(E, "xbar", sm, ordering)
    else
        return 0
    end
end

function Fmass_tint_func(i, sm, ordering, nubar)
    es = range(0, 100, 2000)
    Fmassint_ne = Interpolations.interpolate((vec(es),), Fmass_old.(es, i, sm, ordering, nubar), Gridded(Linear()))
    return extrapolate(Fmassint_ne, 0.0)
end

Fm_1_112sm_NO_nu = Fmass_tint_func(1, "small", "NO", false)
Fm_1_112sm_NO_nubar = Fmass_tint_func(1, "small", "NO", true)
Fm_1_112sm_IO_nu = Fmass_tint_func(1, "small", "IO", false)
Fm_1_112sm_IO_nubar = Fmass_tint_func(1, "small", "IO", true)
Fm_1_270sm_NO_nu = Fmass_tint_func(1, "large", "NO", false)
Fm_1_270sm_NO_nubar = Fmass_tint_func(1, "large", "NO", true)
Fm_1_270sm_IO_nu = Fmass_tint_func(1, "large", "IO", false)
Fm_1_270sm_IO_nubar = Fmass_tint_func(1, "large", "IO", true)
Fm_1_bh_NO_nu = Fmass_tint_func(1, "bh", "NO", false)
Fm_1_bh_NO_nubar = Fmass_tint_func(1, "bh", "NO", true)
Fm_1_bh_IO_nu = Fmass_tint_func(1, "bh", "IO", false)
Fm_1_bh_IO_nubar = Fmass_tint_func(1, "bh", "IO", true)

Fm_2_112sm_NO_nu = Fmass_tint_func(2, "small", "NO", false)
Fm_2_112sm_NO_nubar = Fmass_tint_func(2, "small", "NO", true)
Fm_2_112sm_IO_nu = Fmass_tint_func(2, "small", "IO", false)
Fm_2_112sm_IO_nubar = Fmass_tint_func(2, "small", "IO", true)
Fm_2_270sm_NO_nu = Fmass_tint_func(2, "large", "NO", false)
Fm_2_270sm_NO_nubar = Fmass_tint_func(2, "large", "NO", true)
Fm_2_270sm_IO_nu = Fmass_tint_func(2, "large", "IO", false)
Fm_2_270sm_IO_nubar = Fmass_tint_func(2, "large", "IO", true)
Fm_2_bh_NO_nu = Fmass_tint_func(2, "bh", "NO", false)
Fm_2_bh_NO_nubar = Fmass_tint_func(2, "bh", "NO", true)
Fm_2_bh_IO_nu = Fmass_tint_func(2, "bh", "IO", false)
Fm_2_bh_IO_nubar = Fmass_tint_func(2, "bh", "IO", true)

Fm_3_112sm_NO_nu = Fmass_tint_func(3, "small", "NO", false)
Fm_3_112sm_NO_nubar = Fmass_tint_func(3, "small", "NO", true)
Fm_3_112sm_IO_nu = Fmass_tint_func(3, "small", "IO", false)
Fm_3_112sm_IO_nubar = Fmass_tint_func(3, "small", "IO", true)
Fm_3_270sm_NO_nu = Fmass_tint_func(3, "large", "NO", false)
Fm_3_270sm_NO_nubar = Fmass_tint_func(3, "large", "NO", true)
Fm_3_270sm_IO_nu = Fmass_tint_func(3, "large", "IO", false)
Fm_3_270sm_IO_nubar = Fmass_tint_func(3, "large", "IO", true)
Fm_3_bh_NO_nu = Fmass_tint_func(3, "bh", "NO", false)
Fm_3_bh_NO_nubar = Fmass_tint_func(3, "bh", "NO", true)
Fm_3_bh_IO_nu = Fmass_tint_func(3, "bh", "IO", false)
Fm_3_bh_IO_nubar = Fmass_tint_func(3, "bh", "IO", true)

function Fmass_test1(E, sm, ordering, nubar)
    if sm == "small"
        if ordering == "NO"
            if nubar == false
                return Fm_1_112sm_NO_nu(E)
            elseif nubar == true
                return Fm_1_112sm_NO_nubar(E)
            else
                return 0.0
            end
        elseif ordering == "IO"
            if nubar == false
                return Fm_1_112sm_IO_nu(E)
            elseif nubar == true
                return Fm_1_112sm_IO_nubar(E)
            else
                return 0.0
            end
        else
            return 0.0
        end
    elseif sm == "large"
        if ordering == "NO"
            if nubar == false
                return Fm_1_270sm_NO_nu(E)
            elseif nubar == true
                return Fm_1_270sm_NO_nubar(E)
            else
                return 0.0
            end
        elseif ordering == "IO"
            if nubar == false
                return Fm_1_270sm_IO_nu(E)
            elseif nubar == true
                return Fm_1_270sm_IO_nubar(E)
            else
                return 0.0
            end
        else
            return 0.0
        end
    elseif sm == "bh"
        if ordering == "NO"
            if nubar == false
                return Fm_1_bh_NO_nu(E)
            elseif nubar == true
                return Fm_1_bh_NO_nubar(E)
            else
                return 0.0
            end
        elseif ordering == "IO"
            if nubar == false
                return Fm_1_bh_IO_nu(E)
            elseif nubar == true
                return Fm_1_bh_IO_nubar(E)
            else
                return 0.0
            end
        else
            return 0.0
        end
    else
        return 0.0
    end
end
function Fmass_test2(E, sm, ordering, nubar)
    if sm == "small"
        if ordering == "NO"
            if nubar == false
                return Fm_2_112sm_NO_nu(E)
            elseif nubar == true
                return Fm_2_112sm_NO_nubar(E)
            else
                return 0.0
            end
        elseif ordering == "IO"
            if nubar == false
                return Fm_2_112sm_IO_nu(E)
            elseif nubar == true
                return Fm_2_112sm_IO_nubar(E)
            else
                return 0.0
            end
        else
            return 0.0
        end
    elseif sm == "large"
        if ordering == "NO"
            if nubar == false
                return Fm_2_270sm_NO_nu(E)
            elseif nubar == true
                return Fm_2_270sm_NO_nubar(E)
            else
                return 0.0
            end
        elseif ordering == "IO"
            if nubar == false
                return Fm_2_270sm_IO_nu(E)
            elseif nubar == true
                return Fm_2_270sm_IO_nubar(E)
            else
                return 0.0
            end
        else
            return 0.0
        end
    elseif sm == "bh"
        if ordering == "NO"
            if nubar == false
                return Fm_2_bh_NO_nu(E)
            elseif nubar == true
                return Fm_2_bh_NO_nubar(E)
            else
                return 0.0
            end
        elseif ordering == "IO"
            if nubar == false
                return Fm_2_bh_IO_nu(E)
            elseif nubar == true
                return Fm_2_bh_IO_nubar(E)
            else
                return 0.0
            end
        else
            return 0.0
        end
    else
        return 0.0
    end
end
function Fmass_test3(E, sm, ordering, nubar)
    if sm == "small"
        if ordering == "NO"
            if nubar == false
                return Fm_3_112sm_NO_nu(E)
            elseif nubar == true
                return Fm_3_112sm_NO_nubar(E)
            else
                return 0.0
            end
        elseif ordering == "IO"
            if nubar == false
                return Fm_3_112sm_IO_nu(E)
            elseif nubar == true
                return Fm_3_112sm_IO_nubar(E)
            else
                return 0.0
            end
        else
            return 0.0
        end
    elseif sm == "large"
        if ordering == "NO"
            if nubar == false
                return Fm_3_270sm_NO_nu(E)
            elseif nubar == true
                return Fm_3_270sm_NO_nubar(E)
            else
                return 0.0
            end
        elseif ordering == "IO"
            if nubar == false
                return Fm_3_270sm_IO_nu(E)
            elseif nubar == true
                return Fm_3_270sm_IO_nubar(E)
            else
                return 0.0
            end
        else
            return 0.0
        end
    elseif sm == "bh"
        if ordering == "NO"
            if nubar == false
                return Fm_3_bh_NO_nu(E)
            elseif nubar == true
                return Fm_3_bh_NO_nubar(E)
            else
                return 0.0
            end
        elseif ordering == "IO"
            if nubar == false
                return Fm_3_bh_IO_nu(E)
            elseif nubar == true
                return Fm_3_bh_IO_nubar(E)
            else
                return 0.0
            end
        else
            return 0.0
        end
    else
        return 0.0
    end
end

function Fmass(E, i, sm, ordering, nubar)
    if i == 1
        return Fmass_test1(E, sm, ordering, nubar)
    elseif i == 2
        return Fmass_test2(E, sm, ordering, nubar)
    elseif i == 3
        return Fmass_test3(E, sm, ordering, nubar)
    else
        return 0.0
    end
end