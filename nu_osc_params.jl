# Using the following work: P. F. de Salas, D. V. Forero, S. Gariazzo, P. Mart ́ınez-Mirav ́e, O. Mena, C. A. Ternes, M. T ́ortola, and J. W. F. Valle. 2020 global reassessment of the neutrino oscillation picture. JHEP, 02:071, 2021

θ12_NO = 34.3*(π/180)
θ23_NO = 49.26*(π/180)
θ13_NO = 8.53*(π/180)
δcp_NO = 194*(π/180)

θ12_IO = 34.3*(π/180)
θ23_IO = 49.46*(π/180)
θ13_IO = 8.58*(π/180)
δcp_IO = 284*(π/180)

U1_NO = [1 0 0;
      0 cos(θ23_NO) sin(θ23_NO);
      0 -sin(θ23_NO) cos(θ23_NO)]

U2_NO = [cos(θ13_NO) 0 sin(θ13_NO)*exp(-im*δcp_NO);
      0 1 0;
      -sin(θ13_NO)*exp(im*δcp_NO) 0 cos(θ13_NO)]

U3_NO = [cos(θ12_NO) sin(θ12_NO) 0;
      -sin(θ12_NO) cos(θ12_NO) 0;
      0 0 1]

U1_IO = [1 0 0;
      0 cos(θ23_IO) sin(θ23_IO);
      0 -sin(θ23_IO) cos(θ23_IO)]

U2_IO = [cos(θ13_IO) 0 sin(θ13_IO)*exp(-im*δcp_IO);
      0 1 0;
      -sin(θ13_IO)*exp(im*δcp_IO) 0 cos(θ13_IO)]

U3_IO = [cos(θ12_IO) sin(θ12_IO) 0;
      -sin(θ12_IO) cos(θ12_IO) 0;
      0 0 1]

U_NO = U1_NO*U2_NO*U3_NO

U_IO = U1_IO*U2_IO*U3_IO

Usqred_NO = real(U_NO .* conj(U_NO))
Usqred_IO = real(U_IO .* conj(U_IO))

function Usqred(ordering)
    if ordering == "NO"
        return Usqred_NO
    elseif ordering == "IO"
        return Usqred_IO
    else
        return 0
    end
end

δm2_21 = 7.5e-5 # eV²
δm2_31_NO = 2.55e-3
δm2_31_IO = -2.45e-3

# m2_1_NO (m2_3_NO) set the absolute mass scale, can't be less than 0
m2_1_NO = 0
m2_2_NO = δm2_21
m2_3_NO = δm2_21/2 + δm2_31_NO

m2_3_IO = 0
m2_1_IO = -δm2_31_IO - δm2_21
m2_2_IO = -δm2_31_IO + δm2_21
