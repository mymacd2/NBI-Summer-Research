# module dsnbIOmodule

# export DSNB_vdecay_3ν_νe_IO_1pc

include("F_calculations.jl")

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


# Effective length function
function leff(x) 
    if x == 0.0
        return 0.0
    else
        (3.844126829887412e-6 + 0.9994711237154343*x - 1.2171757160998375*x^2 + 1.2011732904692138*x^3 - 1.0258795580085436*x^4 + 0.7535135022231549*x^5 - 0.4576091931617817*x^6 
        + 0.2207326320476638*x^7 - 0.0820545498767346*x^8 + 0.022953521018636694*x^9 - 0.004712524469463747*x^(10) + 0.0006858791896371295*x^(11) - 6.682066367256103e-5*x^(12) + 3.9015397017829775e-6*x^(13) - 1.0310563080271216e-7*x^(14))/70
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
end


# Invisible decay implementation

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


# Energy spectrum function for SH case
function ψSH(Eh, El, hc)
    if hc
        return 2*El/Eh^2
    else
        return (2/Eh)*(1-(El/Eh))
    end
end


# 3ν IO treatment (2 → 1 QD, 2/1 → 3 SH)

# Case A
# ν2 → ν1, ν2 → ν3, no ν1 → ν3:
# B_21 = 0.5 (no hf), B_23 = 0.25, B_13 = 0
# NB: here if we want to only consider visible decays, we should set α1 to 0 always

# Case B
# ν2 → ν1, ν1 → ν3, no ν2 → ν3:
# B_21 = 1 (no hf), B_13 = 0.5, B_23 = 0

# Case C
# ν2 → ν3, ν1 → ν3, no ν2 → ν1
# B_23 = 0.5, B_13 = 0.5, B_21 = 0

# Democratic
# ν2 → ν3, ν1 → ν3, ν2 → ν1
# B_21 = B_23 = 1/3, B_13 = 0.5

# Note: supposing validity of SH and QD approximations, branching ratios are set in B and C

function branching_3ν_IO(casechoice, j, jbar, i, ibar)
    if casechoice == "A"
        if j == 2 && i == 1 && jbar == ibar
            return 0.5
        elseif j == 2 && i == 3
            return 0.25
        else
            return 0.0
        end
    elseif casechoice == "B"
        if j == 2 && i == 1 && jbar == ibar
            return 1.0
        elseif j == 1 && i == 3
            return 0.5
        else 
            return 0.0
        end
    elseif casechoice == "C"
        if j == 2 && i == 3
            return 0.5
        elseif j == 1 && i == 3
            return 0.5
        else
            return 0.0
        end
    elseif casechoice == "democratic"
        if (j == 2 && i == 1 && jbar == ibar) || (j == 2 && i == 3)
            return 1/3
        elseif j == 1 && i == 3
            return 0.5
        else
            return 0.0
        end
    else
        return 0.0
    end
end

function q21contrib_IO_1pc(E, z0, z, α2, nubar, casechoice, sm, normchoice)

    Ers = E*(1+z)/(1+z0)

    qnorm = 3.086e19 * 1.516e15 * 1e6 / (3e8 * 1e12)
    return qnorm * (c0/Hubble(z)) *  DSNB_idecay_1pc(Ers, z, α2, 2, nubar, "IO", sm, normchoice) * (α2 * branching_3ν_IO(casechoice, 2, nubar, 1, nubar)/Ers)
end

function DSNB_vdecay_1_IO_1pc(E, z0, α2, α1, nubar, casechoice, sm, normchoice)
    if casechoice == "A"
        α1 = 0
    end
    integrand(z) = (DSNB_integrand_1pc(E, z0, z, 1, nubar, "IO", sm, normchoice) +
                        q21contrib_IO_1pc(E, z0, z, α2, nubar, casechoice, sm, normchoice))*decay(E, α1, z0, z)
    zcut = zcutoff(E, α1, z0)
    if zcut > z0 && zcut < 5 && decay(E, α1, z0, zcut) < 0.1
        return (1/(1+z0)) * quadgk(z -> integrand(z), z0, zcut, 5, rtol=1e-2)[1]
    else
        return (1/(1+z0)) * quadgk(z -> integrand(z), z0, 5, rtol=1e-2)[1]
    end
end

function q23contrib_IO_1pc(E, z0, z, α2, twobar, threebar, casechoice, sm, normchoice)

    if twobar == threebar
        hc = true
    else
        hc = false
    end

    Ers = E*(1+z)/(1+z0)

    qnorm = 3.086e19 * 1.516e15 * 1e6 / (3e8 * 1e12)
    integrand(Eprime, z) = qnorm * (c0/Hubble(z))*DSNB_idecay_1pc(Eprime, z, α2, 2, twobar, "IO", sm, normchoice) * (α2 * branching_3ν_IO(casechoice, 2, twobar, 3, threebar)/Eprime) * ψSH(Eprime, Ers, hc)

    Emax = Ers + 50

    Ecutoff = 1.5*100^(1/3)*Ers
    if Ecutoff < Emax && (ψSH(Ecutoff, Ers, hc)/Ecutoff)/(ψSH(Ers, Ers, hc)/Ers) < 0.1
        return quadgk(Eprime -> integrand(Eprime, z), Ers, Ecutoff, Emax, rtol=1e-2)[1]
    else
        return quadgk(Eprime -> integrand(Eprime, z), Ers, Emax, rtol=1e-2)[1]
    end
end

function q13contrib_IO_1pc(E, z0, z, α2, α1, onebar, threebar, casechoice, sm, normchoice)

    if onebar == threebar
        hc = true
    else
        hc = false
    end

    Ers = E*(1+z)/(1+z0)

    qnorm = 3.086e19 * 1.516e15 * 1e6 / (3e8 * 1e12)
    integrand(Eprime, z) = qnorm * (c0/Hubble(z))*DSNB_vdecay_1_IO_1pc(Eprime, z, α2, α1, onebar, casechoice, sm, normchoice) * (α1 * branching_3ν_IO(casechoice, 1, onebar, 3, threebar)/Eprime) * ψSH(Eprime, Ers, hc)

    Emax = Ers + 50

    Ecutoff = 1.5*100^(1/3)*Ers
    if Ecutoff < Emax && (ψSH(Ecutoff, Ers, hc)/Ecutoff)/(ψSH(Ers, Ers, hc)/Ers) < 0.1
        return quadgk(Eprime -> integrand(Eprime, z), Ers, Ecutoff, Emax, rtol=1e-2)[1]
    else
        return quadgk(Eprime -> integrand(Eprime, z), Ers, Emax, rtol=1e-2)[1]
    end
end

function DSNB_vdecay_3_IO_1pc(E, α2, α1, threebar, casechoice, sm, normchoice)
    integrand(z) = (DSNB_integrand_1pc(E, 0, z, 3, threebar, "IO", sm, normchoice)
                    + q23contrib_IO_1pc(E, 0, z, α2, true, threebar, casechoice, sm, normchoice)
                    + q23contrib_IO_1pc(E, 0, z, α2, false, threebar, casechoice, sm, normchoice)
                    + q13contrib_IO_1pc(E, 0, z, α2, α1, true, threebar, casechoice, sm, normchoice)
                    + q13contrib_IO_1pc(E, 0, z, α2, α1, false, threebar, casechoice, sm, normchoice))
    return quadgk(z -> integrand(z), 0, 5, rtol=1e-2)[1]
end

function DSNB_vdecay_3ν_νe_IO_1pc(E, α2, α1, nubar, casechoice, sm, normchoice)

    ν3 = DSNB_vdecay_3_IO_1pc(E, α2, α1, nubar, casechoice, sm, normchoice)
    ν2 = DSNB_idecay_1pc(E, 0, α2, 2, nubar, "IO", sm, normchoice)
    ν1 = DSNB_vdecay_1_IO_1pc(E, 0, α2, α1, nubar, casechoice, sm, normchoice)
    
    return Usqred("IO")[1, 1]*ν1 + Usqred("IO")[1, 2]*ν2 + Usqred("IO")[1, 3]*ν3
end;

# end