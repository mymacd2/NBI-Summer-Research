using DelimitedFiles
using Interpolations
using QuadGK
using SpecialFunctions

c0 = 3e8 # m s^(-1)
H0 = 70.0 # km s^(-1) Mpc^(-1)

# Star formation rate (SFR)
SFR(z) = ((1 + z)^(-34) + ((1 + z)/5000)^3 + ((1 + z)/9)^(35))^(-0.1)

#Initial mass function (IMF)
η(M) = M^(-2.35)

# Unnormalized supernova rate (SNR0)
SNR0(z, M) = η(M) * SFR(z)

intSNR0 = quadgk(M -> SNR0(0, M), 8, 125)[1]
SNRnorm = 1.25e-4 / intSNR0 # units: 1.25e-4 Mpc^(-3) yr^(-1)

# Normalized SFR
SNR(z, M) = SNRnorm*SNR0(z, M);

# Hubble parameter dep on redshift:
function Hubble(z) 
    H0 = 70 # km s^(-1) Mpc^(-1)
    energy_matter, energy_dark = 0.3, 0.7
    return H0*sqrt(energy_matter*((1+z)^3) + energy_dark)
end

# Antiderivative of IMF evaluated at endpts
ηAD(m_max, m_min) = (m_max^(-1.35)-m_min^(-1.35))/(-1.35)