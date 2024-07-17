# Running this file will update/add the file "F0s_vec.txt", which will have initial fluxes evaluated on an energy grid that is (currently) running from 0 to 100 MeV. Simply interpolate on the grid

using DelimitedFiles
using Plots
using Interpolations
using LaTeXStrings
using QuadGK
using SpecialFunctions

#              time        luminosity               <e>             <e^2>
#               "s"           "foe/s"           "MeV^1"           "MeV^2"

νe_112sm = readdlm("Data/neutrino_signal_nu_e-LS220-s11.2c.data", comments=true)
νebar_112sm = readdlm("Data/neutrino_signal_nubar_e-LS220-s11.2c.data", comments=true)
νx_112sm = readdlm("Data/neutrino_signal_nu_x-LS220-s11.2c.data", comments=true)

νe_270sm = readdlm("Data/neutrino_signal_nu_e-s27.0c-LS220.data", comments=true)
νebar_270sm = readdlm("Data/neutrino_signal_nubar_e-s27.0c-LS220.data", comments=true)
νx_270sm = readdlm("Data/neutrino_signal_nu_x-s27.0c-LS220.data", comments=true)

νe_bh = readdlm("Data/s40s7b2c_neutrino_signal_nu_e.dat", comments=true)
νebar_bh = readdlm("Data/s40s7b2c_neutrino_signal_nubar_e.dat", comments=true)
νx_bh = readdlm("Data/s40s7b2c_neutrino_signal_nu_x.dat", comments=true)


# 27 Solar Mass
lum_νe_270sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νe_270sm[:,1]),), νe_270sm[:,2], Gridded(Linear()))
lum_νe_270sm = extrapolate(lum_νe_270sm_ne, 0.0)

em1_νe_270sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νe_270sm[:,1]),), νe_270sm[:,3], Gridded(Linear()))
em1_νe_270sm = extrapolate(em1_νe_270sm_ne, 0.0)

em2_νe_270sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νe_270sm[:,1]),), νe_270sm[:,4], Gridded(Linear()))
em2_νe_270sm = extrapolate(em2_νe_270sm_ne, 0.0)

lum_νebar_270sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νebar_270sm[:,1]),), νebar_270sm[:,2], Gridded(Linear()))
lum_νebar_270sm = extrapolate(lum_νebar_270sm_ne, 0.0)

em1_νebar_270sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νebar_270sm[:,1]),), νebar_270sm[:,3], Gridded(Linear()))
em1_νebar_270sm = extrapolate(em1_νebar_270sm_ne, 0.0)

em2_νebar_270sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νebar_270sm[:,1]),), νebar_270sm[:,4], Gridded(Linear()))
em2_νebar_270sm = extrapolate(em2_νebar_270sm_ne, 0.0)

lum_νx_270sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νx_270sm[:,1]),), νx_270sm[:,2], Gridded(Linear()))
lum_νx_270sm = extrapolate(lum_νx_270sm_ne, 0.0)

em1_νx_270sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νx_270sm[:,1]),), νx_270sm[:,3], Gridded(Linear()))
em1_νx_270sm = extrapolate(em1_νx_270sm_ne, 0.0)

em2_νx_270sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νx_270sm[:,1]),), νx_270sm[:,4], Gridded(Linear()))
em2_νx_270sm = extrapolate(em2_νx_270sm_ne, 0.0)

# 11.2 Solar Mass
lum_νe_112sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νe_112sm[:,1]),), νe_112sm[:,2], Gridded(Linear()))
lum_νe_112sm = extrapolate(lum_νe_112sm_ne, 0.0)

em1_νe_112sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νe_112sm[:,1]),), νe_112sm[:,3], Gridded(Linear()))
em1_νe_112sm = extrapolate(em1_νe_112sm_ne, 0.0)

em2_νe_112sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νe_112sm[:,1]),), νe_112sm[:,4], Gridded(Linear()))
em2_νe_112sm = extrapolate(em2_νe_112sm_ne, 0.0)

lum_νebar_112sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νebar_112sm[:,1]),), νebar_112sm[:,2], Gridded(Linear()))
lum_νebar_112sm = extrapolate(lum_νebar_112sm_ne, 0.0)

em1_νebar_112sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νebar_112sm[:,1]),), νebar_112sm[:,3], Gridded(Linear()))
em1_νebar_112sm = extrapolate(em1_νebar_112sm_ne, 0.0)

em2_νebar_112sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νebar_112sm[:,1]),), νebar_112sm[:,4], Gridded(Linear()))
em2_νebar_112sm = extrapolate(em2_νebar_112sm_ne, 0.0)

lum_νx_112sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νx_112sm[:,1]),), νx_112sm[:,2], Gridded(Linear()))
lum_νx_112sm = extrapolate(lum_νx_112sm_ne, 0.0)

em1_νx_112sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νx_112sm[:,1]),), νx_112sm[:,3], Gridded(Linear()))
em1_νx_112sm = extrapolate(em1_νx_112sm_ne, 0.0)

em2_νx_112sm_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νx_112sm[:,1]),), νx_112sm[:,4], Gridded(Linear()))
em2_νx_112sm = extrapolate(em2_νx_112sm_ne, 0.0)

# Black Hole
lum_νe_bh_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νe_bh[:,1]),), νe_bh[:,2], Gridded(Linear()))
lum_νe_bh = extrapolate(lum_νe_bh_ne, 0.0)

em1_νe_bh_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νe_bh[:,1]),), νe_bh[:,3], Gridded(Linear()))
em1_νe_bh = extrapolate(em1_νe_bh_ne, 0.0)

em2_νe_bh_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νe_bh[:,1]),), νe_bh[:,4], Gridded(Linear()))
em2_νe_bh = extrapolate(em2_νe_bh_ne, 0.0)

lum_νebar_bh_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νebar_bh[:,1]),), νebar_bh[:,2], Gridded(Linear()))
lum_νebar_bh = extrapolate(lum_νebar_bh_ne, 0.0)

em1_νebar_bh_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νebar_bh[:,1]),), νebar_bh[:,3], Gridded(Linear()))
em1_νebar_bh = extrapolate(em1_νebar_bh_ne, 0.0)

em2_νebar_bh_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νebar_bh[:,1]),), νebar_bh[:,4], Gridded(Linear()))
em2_νebar_bh = extrapolate(em2_νebar_bh_ne, 0.0)

lum_νx_bh_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νx_bh[:,1]),), νx_bh[:,2], Gridded(Linear()))
lum_νx_bh = extrapolate(lum_νx_bh_ne, 0.0)

em1_νx_bh_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νx_bh[:,1]),), νx_bh[:,3], Gridded(Linear()))
em1_νx_bh = extrapolate(em1_νx_bh_ne, 0.0)

em2_νx_bh_ne = Interpolations.interpolate((Interpolations.deduplicate_knots!(νx_bh[:,1]),), νx_bh[:,4], Gridded(Linear()))
em2_νx_bh = extrapolate(em2_νx_bh_ne, 0.0)



# β = e (electron), ebar (anti electron), x (non electron)
# sm = small (11.2 SM), large (27.0 SM), bh (Black Hole)

function L(t, β, sm)
    if β == "e" && sm == "small"
        return lum_νe_112sm(t)
    elseif β == "e" && sm == "large"
        return lum_νe_270sm(t)
    elseif β == "e" && sm == "bh"
        return lum_νe_bh(t)
    elseif β == "ebar" && sm == "small"
        return lum_νebar_112sm(t)
    elseif β == "ebar" && sm == "large"
        return lum_νebar_270sm(t)
    elseif β == "ebar" && sm == "bh"
        return lum_νebar_bh(t)    
    elseif β == "x" && sm == "small"
        return lum_νx_112sm(t)
    elseif β == "x" && sm == "large"
        return lum_νx_270sm(t)
    elseif β == "x" && sm == "bh"
        return lum_νx_bh(t)
    else
        return 0
    end
end

function Em1(t, β, sm)
    if β == "e" && sm == "small"
        return em1_νe_112sm(t)
    elseif β == "e" && sm == "large"
        return em1_νe_270sm(t)
    elseif β == "e" && sm == "bh"
        return em1_νe_bh(t)
    elseif β == "ebar" && sm == "small"
        return em1_νebar_112sm(t)
    elseif β == "ebar" && sm == "large"
        return em1_νebar_270sm(t)
    elseif β == "ebar" && sm == "bh"
        return em1_νebar_bh(t)    
    elseif β == "x" && sm == "small"
        return em1_νx_112sm(t)
    elseif β == "x" && sm == "large"
        return em1_νx_270sm(t)
    elseif β == "x" && sm == "bh"
        return em1_νx_bh(t)
    else
        return 0
    end
end

function Em2(t, β, sm)
    if β == "e" && sm == "small"
        return em2_νe_112sm(t)
    elseif β == "e" && sm == "large"
        return em2_νe_270sm(t)
    elseif β == "e" && sm == "bh"
        return em2_νe_bh(t)
    elseif β == "ebar" && sm == "small"
        return em2_νebar_112sm(t)
    elseif β == "ebar" && sm == "large"
        return em2_νebar_270sm(t)
    elseif β == "ebar" && sm == "bh"
        return em2_νebar_bh(t)    
    elseif β == "x" && sm == "small"
        return em2_νx_112sm(t)
    elseif β == "x" && sm == "large"
        return em2_νx_270sm(t)
    elseif β == "x" && sm == "bh"
        return em2_νx_bh(t)
    else
        return 0
    end
end


# Pinching parameter expression

α(t, β, sm) = (2*Em1(t, β, sm)^2 - Em2(t, β, sm))/(Em2(t, β, sm) - Em1(t, β, sm)^2)


# (Normalized) energy distribution

# Unnormalized first
function ϕ0(E, t, β, sm) 

    if Em1(t, β, sm) == 0
        return 0.0
    else
        return (E/Em1(t, β, sm))^(α(t, β, sm))*exp(-(E*(1+α(t, β, sm)))/Em1(t, β, sm))
    end
end

intϕ_exact(t, β, sm) = Em1(t, β, sm)^(-α(t, β, sm))*((1+α(t, β, sm))/Em1(t, β, sm))^(-1-α(t, β, sm))*gamma(1+α(t, β, sm))

# Using analytical normalizaiton
function ϕ(E, t, β, sm)
    if Em1(t, β, sm) == 0
        return 0.0
    else
        return ϕ0(E, t, β, sm) / intϕ_exact(t, β, sm)
    end
end


# Calculating the initial (unoscillated) SN fluxes

function dF0(E, t, β, sm) 
    if Em1(t, β, sm) == 0
        return 0.0
    else
        return L(t, β, sm) * ϕ(E, t, β, sm)/Em1(t, β, sm)
    end
end

loe_to_MeV = 6.2415e56

function fake_tint_F0(energy, β, sm)
    ts_int = range(-0.4, 9, 10000)
    dts_int = ts_int[2] - ts_int[1]
    return sum(dts_int * dF0.(energy, ts_int, β, sm))
end

function F0_tint(β, sm)
    es = range(0, 50, 1000)
    tintF0 = fake_tint_F0.(es, β, sm)
    F0int_ne = Interpolations.interpolate((vec(es),), tintF0, Gridded(Linear()))
    return extrapolate(F0int_ne, 0.0)
end

function F0_tint_vec(β, sm)
    es = range(0, 100, 2000)
    tintF0 = fake_tint_F0.(es, β, sm)
    return tintF0
end

# F0_νe_270sm_vec, F0_νebar_270sm_vec, F0_νx_270sm_vec = F0_tint_vec("e", "large"), F0_tint_vec("ebar", "large"), F0_tint_vec("x", "large")
# F0_νe_112sm_vec, F0_νebar_112sm_vec, F0_νx_112sm_vec = F0_tint_vec("e", "small"), F0_tint_vec("ebar", "small"), F0_tint_vec("x", "small")
# F0_νe_bh_vec, F0_νebar_bh_vec, F0_νx_bh_vec = F0_tint_vec("e", "bh"), F0_tint_vec("ebar", "bh"), F0_tint_vec("x", "bh");

# F0s_vec = hcat(F0_νe_270sm_vec, F0_νebar_270sm_vec, F0_νx_270sm_vec, F0_νe_112sm_vec, F0_νebar_112sm_vec, F0_νx_112sm_vec, F0_νe_bh_vec, F0_νebar_bh_vec, F0_νx_bh_vec)

# writedlm("F0s_vec.txt", F0s_vec)