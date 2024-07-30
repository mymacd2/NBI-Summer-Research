include("F_calculations.jl")

# We can take normchoice to be SNRnorm, SNRnorm_low, or SNRnorm_high
function DSNB_integrand(E, z0, z, i, nubar, ordering, bh_frac, normchoice)
    c0 = 3e8 # m s^(-1)
    if bh_frac == "21"
        return 0.00208 * (c0*normchoice*SFR(z)/Hubble(z)) * (ηAD(15, 8)*Fmass(E*(1+z)/(1+z0), i, "small", ordering, nubar) 
        + (ηAD(22, 15)+ηAD(27, 25))*Fmass(E*(1+z)/(1+z0), i, "large", ordering, nubar) + (ηAD(25, 22)+ηAD(125, 27))*Fmass(E*(1+z)/(1+z0), i, "bh", ordering, nubar))
    elseif bh_frac == "41"
        return 0.00208 * (c0*normchoice*SFR(z)/Hubble(z)) * (ηAD(15, 8)*Fmass(E*(1+z)/(1+z0), i, "small", ordering, nubar) 
        + ηAD(125, 15)*Fmass(E*(1+z)/(1+z0), i, "bh", ordering, nubar))
    elseif bh_frac == "09"
        return 0.00208 * (c0*normchoice*SFR(z)/Hubble(z)) * (ηAD(15, 8)*Fmass(E*(1+z)/(1+z0), i, "small", ordering, nubar)
        + ηAD(40, 15)*Fmass(E*(1+z)/(1+z0), i, "large", ordering, nubar) + ηAD(125, 40)*Fmass(E*(1+z)/(1+z0), i, "bh", ordering, nubar))
    else
        return 0
    end
end

# Mass method
DSNB(E, z0, i, nubar, ordering, bh_frac, normchoice) = 1/(1+z0) * quadgk(z -> DSNB_integrand(E, z0, z, i, nubar, ordering, bh_frac, normchoice), z0, 5, rtol=1e-2)[1]

# νe method
function DSNB(E, z0, nubar, ordering, bh_frac, normchoice)
    ν1 = DSNB(E, z0, 1, nubar, ordering, bh_frac, normchoice)
    ν2 = DSNB(E, z0, 2, nubar, ordering, bh_frac, normchoice)
    ν3 = DSNB(E, z0, 3, nubar, ordering, bh_frac, normchoice)
    return Usqred(ordering)[1, 1]*ν1 + Usqred(ordering)[1, 2]*ν2 + Usqred(ordering)[1, 3]*ν3
end


# Separating the progenitors: 1pc ≡ one progenitor contribution

function DSNB_integrand_1pc(E, z0, z, i, nubar, ordering, sm, normchoice)
    c0 = 3e8 # m s^(-1)
    # Adding in this integrated mass function so we can multiply by bh fractions at the end
    return 0.00208 * ηAD(125, 8) * (c0*normchoice*SFR(z)/Hubble(z)) * Fmass(E*(1+z)/(1+z0), i, sm, ordering, nubar)
end

DSNB_1pc(E, z0, i, nubar, ordering, sm, normchoice) = 1/(1+z0) * quadgk(z -> DSNB_integrand_1pc(E, z0, z, i, nubar, ordering, sm, normchoice), z0, 5, rtol=1e-2)[1]

function DSNB_1pc(E, z0, nubar, ordering, sm, normchoice)
    ν1 = DSNB_1pc(E, z0, 1, nubar, ordering, sm, normchoice)
    ν2 = DSNB_1pc(E, z0, 2, nubar, ordering, sm, normchoice)
    ν3 = DSNB_1pc(E, z0, 3, nubar, ordering, sm, normchoice)
    return Usqred(ordering)[1, 1]*ν1 + Usqred(ordering)[1, 2]*ν2 + Usqred(ordering)[1, 3]*ν3
end;


# Define this effective length function
LeffIntegrand(z) = (H0*sqrt(energy_matter*(1+z)^3 + energy_dark))^(-1)*(1+z)^(-2)

leff_quadint(z0) = quadgk(z -> LeffIntegrand(z), 0, z0)[1]

zs_leff = range(0, 5, 1000)
leff_ne = Interpolations.interpolate((vec(zs_leff),), leff_quadint.(zs_leff), Gridded(Linear()))
lefftrue = extrapolate(leff_ne, 0.0)

# Approximating with a 14 degree polynomial

# Procedure for fitting the polynomial
#=
xsvec = range(0, 5, 100)
leffsvec = 70*leff.(xsvec)

p = Polynomials.fit(xsvec, leffsvec, 14)
=#

function leff(x) 
    if x == 0.0
        return 0.0
    else
        (3.844126829887412e-6 + 0.9994711237154343*x - 1.2171757160998375*x^2 + 1.2011732904692138*x^3 - 1.0258795580085436*x^4 + 0.7535135022231549*x^5 - 0.4576091931617817*x^6 + 0.2207326320476638*x^7 - 0.0820545498767346*x^8 + 0.022953521018636694*x^9 - 0.004712524469463747*x^(10) + 0.0006858791896371295*x^(11) - 6.682066367256103e-5*x^(12) + 3.9015397017829775e-6*x^(13) - 1.0310563080271216e-7*x^(14))/70
    end
end

# Decay function 
function decay(E, α, z0, z)
    scalefactor = 4.68e28
    return exp(-scalefactor*α*(leff(z) - leff(z0))*(1+z0)/E)
end


# This function returns the redshift point where the decay function has attenuated by 
# roughly two orders of magnitude, and if we break the integral bounds at this point the integral might run faster
# with QuadGK.jl
function zcutoff(E, α, z0)
    innerpart = 0.007+(E/(4.68e28*α))*log(0.01*exp(-4.68e28*α*leff(z0)/E))
    if innerpart > 0.0
        return -0.555*(log(350*(innerpart))-0.9)
    else
        return 0.0
    end
end;


# Invisible decay implementation

# Mass method
function DSNB_idecay(E, z0, α, i, nubar, ordering, bh_frac, normchoice)
    zcut = zcutoff(E, α, z0)
    if zcut > z0 && zcut < 5 && decay(E, α, z0, zcut) < 0.1
        return quadgk(z -> 1/(1+z0) * DSNB_integrand(E, z0, z, i, nubar, ordering, bh_frac, normchoice)*decay(E, α, z0, z), z0, zcut, 5, rtol=1e-2)[1]
    else
        return quadgk(z -> 1/(1+z0) * DSNB_integrand(E, z0, z, i, nubar, ordering, bh_frac, normchoice)*decay(E, α, z0, z), z0, 5, rtol=1e-2)[1]
    end
end

# νe method
function DSNB_idecay(E, z0, α1, α2, α3, nubar, ordering, bh_frac, normchoice)
    ν1 = DSNB_idecay(E, z0, α1, 1, nubar, ordering, bh_frac, normchoice)
    ν2 = DSNB_idecay(E, z0, α2, 2, nubar, ordering, bh_frac, normchoice)
    ν3 = DSNB_idecay(E, z0, α3, 3, nubar, ordering, bh_frac, normchoice)
    return Usqred(ordering)[1, 1]*ν1 + Usqred(ordering)[1, 2]*ν2 + Usqred(ordering)[1, 3]*ν3
end

# Single progenitor versions
function DSNB_idecay_1pc(E, z0, α, i, nubar, ordering, sm, normchoice)
    zcut = zcutoff(E, α, z0)
    if zcut > z0 && zcut < 5 && decay(E, α, z0, zcut) < 0.1
        return quadgk(z -> 1/(1+z0) * DSNB_integrand_1pc(E, z0, z, i, nubar, ordering, sm, normchoice)*decay(E, α, z0, z), z0, zcut, 5, rtol=1e-2)[1]
    else
        return quadgk(z -> 1/(1+z0) * DSNB_integrand_1pc(E, z0, z, i, nubar, ordering, sm, normchoice)*decay(E, α, z0, z), z0, 5, rtol=1e-2)[1]
    end
end

function DSNB_idecay_1pc(E, z0, α1, α2, α3, nubar, ordering, sm, normchoice)
    ν1 = DSNB_idecay_1pc(E, z0, α1, 1, nubar, ordering, sm, normchoice)
    ν2 = DSNB_idecay_1pc(E, z0, α2, 2, nubar, ordering, sm, normchoice)
    ν3 = DSNB_idecay_1pc(E, z0, α3, 3, nubar, ordering, sm, normchoice)
    return Usqred(ordering)[1, 1]*ν1 + Usqred(ordering)[1, 2]*ν2 + Usqred(ordering)[1, 3]*ν3
end


# Visible 2ν decay implementation:

function qcontrib_2ν_1pc(E, z0, z, j, jbar, i, ibar, αj, ordering, sm, normchoice)

    if ordering == "NO" && j <= i
        return println("error: not kinematically allowed")
    elseif ordering == "IO" && j == 1 && i == 2
        return println("error: not kinematically allowed")
    elseif ordering == "IO" && j == 3 && i == 1
        return println("error: not kinematically allowed")
    elseif ordering == "IO" && j == 3 && i == 2
        return println("error: not kinematically allowed")
    elseif ordering == "IO" && j == 2 && i == 1
        if jbar != ibar
            return 0.0
        else
            Ers = E*(1+z)/(1+z0)
            qnorm = 3.086e19 * 1.516e15 * 1e6 / (3e8 * 1e12)
            return qnorm * (c0/Hubble(z))*DSNB_idecay_1pc(Ers, z, αj, j, jbar, ordering, sm, normchoice) * (αj/Ers)
        end
    elseif ordering == "NOQD" && j > i
        if jbar != ibar
            return 0.0
        else
            Ers = E*(1+z)/(1+z0)
            qnorm = 3.086e19 * 1.516e15 * 1e6 / (3e8 * 1e12)
            return qnorm * (c0/Hubble(z))*DSNB_idecay_1pc(Ers, z, αj, j, jbar, "NO", sm, normchoice) * (αj/Ers)
        end
    else
        if jbar == ibar
            hc = true
        else
            hc = false
        end

        Ers = E*(1+z)/(1+z0)

        qnorm = 3.086e19 * 1.516e15 * 1e6 / (3e8 * 1e12)
        integrand(Eprime) = qnorm * (c0/Hubble(z))*DSNB_idecay_1pc(Eprime, z, αj, j, jbar, ordering, sm, normchoice) * (αj * 0.5/Eprime) * ψSH(Eprime, Ers, hc)
        Emax = Ers + 50

        return quadgk(Eprime -> integrand(Eprime), Ers, Emax, rtol=1e-2)[1]
    end
end

function DSNB_vdecay_2ν_1pc(E, αj, i, ibar, ordering, sm, normchoice)
    if αj == 0
        if ordering == "NO" || ordering == "NOQD"
            return DSNB_idecay_1pc(E, 0, 0, i, ibar, "NO", sm, normchoice)
        else
            return DSNB_idecay_1pc(E, 0, 0, i, ibar, ordering, sm, normchoice)
        end
    elseif ordering == "NO" && i == 3
        return DSNB_idecay_1pc(E, 0, αj, i, ibar, "NO", sm, normchoice)
    elseif ordering == "NOQD" && i == 3
        return DSNB_idecay_1pc(E, 0, αj, i, ibar, "NO", sm, normchoice)
    elseif ordering == "IO" && i == 2
        return DSNB_idecay_1pc(E, 0, αj, i, ibar, "IO", sm, normchoice)
    else
        # Here, for NO we take ν2 to be stable and for IO we take ν1 to be stable, and we take the lightest mass states to be stable as well
        if ordering == "NO" || ordering == "NOQD"
            qint = quadgk(z -> (qcontrib_2ν_1pc(E, 0, z, 3, false, i, ibar, αj, ordering, sm, normchoice)+qcontrib_2ν_1pc(E, 0, z, 3, true, i, ibar, αj, ordering, sm, normchoice)), 0, 5, rtol=1e-2)[1]
            return DSNB_idecay_1pc(E, 0, 0, i, ibar, "NO", sm, normchoice) + qint
        elseif ordering == "IO"
            qint = quadgk(z -> (qcontrib_2ν_1pc(E, 0, z, 2, false, i, ibar, αj, "IO", sm, normchoice)+qcontrib_2ν_1pc(E, 0, z, 2, true, i, ibar, αj, "IO", sm, normchoice)), 0, 5, rtol=1e-2)[1]
            return DSNB_idecay_1pc(E, 0, 0, i, ibar, ordering, sm, normchoice) + qint
        else
            return println("error: ordering takes either 'NO' or 'IO'")
        end
    end
end;

# Other few cases where the heaviest mass state and lightest mass state are stable, and the middle mass state decays into the lightest mass state
function DSNB_vdecay_2ν_alt_1pc(E, αj, i, ibar, ordering, sm, normchoice)
    if αj == 0
        if ordering == "NO" || ordering == "NOQD"
            return DSNB_idecay_1pc(E, 0, 0, i, ibar, "NO", sm, normchoice)
        else
            return DSNB_idecay_1pc(E, 0, 0, i, ibar, ordering, sm, normchoice)
        end
    elseif ordering == "NO" && i == 2
        return DSNB_idecay_1pc(E, 0, αj, i, ibar, "NO", sm, normchoice)
    elseif ordering == "NOQD" && i == 2
        return DSNB_idecay_1pc(E, 0, αj, i, ibar, "NO", sm, normchoice)
    elseif ordering == "IO" && i == 1
        return DSNB_idecay_1pc(E, 0, αj, i, ibar, "IO", sm, normchoice)
    elseif ordering == "NO" && i == 3
        return DSNB_idecay_1pc(E, 0, 0, i, ibar, "NO", sm, normchoice)
    elseif ordering == "IO" && i == 2
        return DSNB_idecay_1pc(E, 0, 0, i, ibar, "IO", sm, normchoice)
    else
        # Here, for NO we take ν1 to be stable and for IO we take ν2 to be stable, and we take the lightest mass states to be stable as well
        if ordering == "NO" || ordering == "NOQD"
            qint = quadgk(z -> (qcontrib_2ν_1pc(E, 0, z, 2, false, i, ibar, αj, ordering, sm, normchoice)+qcontrib_2ν_1pc(E, 0, z, 2, true, i, ibar, αj, ordering, sm, normchoice)), 0, 5, rtol=1e-2)[1]
            return DSNB_idecay_1pc(E, 0, 0, i, ibar, "NO", sm, normchoice) + qint
        elseif ordering == "IO"
            qint = quadgk(z -> (qcontrib_2ν_1pc(E, 0, z, 1, false, i, ibar, αj, "IO", sm, normchoice)+qcontrib_2ν_1pc(E, 0, z, 1, true, i, ibar, αj, "IO", sm, normchoice)), 0, 5, rtol=1e-2)[1]
            return DSNB_idecay_1pc(E, 0, 0, i, ibar, ordering, sm, normchoice) + qint
        else
            return println("error: ordering takes either 'NO' or 'IO'")
        end
    end
end

# daughter specifies the state the heaviest mass state decays into in our 2ν framework
function DSNB_vdecay_2ν_νe_1pc(E, α, daughter, ebar, ordering, sm, normchoice)
    if ordering == "NO" || ordering == "NOQD"
        if daughter == 1
            α1, α2, α3 = α, 0.0, α
        elseif daughter == 2
            α1, α2, α3 = 0.0, α, α
        else
            return println("error: for NO, 'daughter' must be either ν1 or ν2")
        end
    elseif ordering == "IO"
        if daughter == 3
            α3, α1, α2 = α, 0.0, α
        elseif daughter == 1
            α3, α1, α2 = 0.0, α, α
        else
            return println("error: for IO, 'daughter' must be either ν1 or ν3")
        end
    else
        return println("error: ordering must take either 'NO' or 'IO'")
    end
    ν3 = DSNB_vdecay_2ν_1pc(E, α3, 3, ebar, ordering, sm, normchoice)
    ν2 = DSNB_vdecay_2ν_1pc(E, α2, 2, ebar, ordering, sm, normchoice)
    ν1 = DSNB_vdecay_2ν_1pc(E, α1, 1, ebar, ordering, sm, normchoice)
    return Usqred("NO")[1, 1]*ν1 + Usqred("NO")[1, 2]*ν2 + Usqred("NO")[1, 3]*ν3
end

# If we don't specify a daughter, the decay channel is assumed to be from the second lightest to the lightest state
function DSNB_vdecay_2ν_νe_1pc(E, α, ebar, ordering, sm, normchoice)
    if ordering == "NO" || ordering == "NOQD"
        α1, α2, α3 = α, α, 0.0
    elseif ordering == "IO"
        α1, α2, α3 = α, 0.0, α
    else
        return println("error: ordering must take either 'NO' or 'IO'")
    end
    ν3 = DSNB_vdecay_2ν_alt_1pc(E, α3, 3, ebar, ordering, sm, normchoice)
    ν2 = DSNB_vdecay_2ν_alt_1pc(E, α2, 2, ebar, ordering, sm, normchoice)
    ν1 = DSNB_vdecay_2ν_alt_1pc(E, α1, 1, ebar, ordering, sm, normchoice)
    return Usqred("NO")[1, 1]*ν1 + Usqred("NO")[1, 2]*ν2 + Usqred("NO")[1, 3]*ν3
end


# Free black hole fraction:

# bh_frac could be anywhere from 0 to 0.41
function DSNB_freefbh(E, z0, nubar, ordering, fbh, normchoice)
    if ordering == "NOQD"
        ordering = "NO"
    end
    dsnb_small = DSNB_1pc(E, z0, nubar, ordering, "small", normchoice)
    dsnb_large = DSNB_1pc(E, z0, nubar, ordering, "large", normchoice)
    dsnb_bh = DSNB_1pc(E, z0, nubar, ordering, "bh", normchoice)

    fsmall = ηAD(15, 8)/ηAD(125, 8)
    if 1 - fsmall - fbh < 0
        fbh = 0.41
    elseif fbh < 0.09
        fbh = 0.09
    end
    flarge = 1 - fsmall - fbh

    return fsmall*dsnb_small + fbh*dsnb_bh + flarge*dsnb_large
end

function DSNB_idecay_freefbh(E, z0, α1, α2, α3, nubar, ordering, fbh, normchoice)
    if ordering == "NOQD"
        ordering = "NO"
    end
    dsnb_small = DSNB_idecay_1pc(E, z0, α1, α2, α3, nubar, ordering, "small", normchoice)
    dsnb_large = DSNB_idecay_1pc(E, z0, α1, α2, α3, nubar, ordering, "large", normchoice)
    dsnb_bh = DSNB_idecay_1pc(E, z0, α1, α2, α3, nubar, ordering, "bh", normchoice)

    fsmall = ηAD(15, 8)/ηAD(125, 8)
    if 1 - fsmall - fbh < 0
        fbh = 0.41
    elseif fbh < 0.09
        fbh = 0.09
    end
    flarge = 1 - fsmall - fbh

    return fsmall*dsnb_small + fbh*dsnb_bh + flarge*dsnb_large
end;

