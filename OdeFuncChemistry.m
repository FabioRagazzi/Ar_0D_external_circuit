function [dydt] = OdeFuncChemistry(t,y,k_Bolsig,stoichiometric_matrix,Te_Bolsig,mue_N,T,Ngas,Gap_length,Vext,R,Area)
    n = y(1:end-1);
    v = y(end);
    ne = n(1);
    nAr = n(2);
    nAr2p = n(3);
    nArp = n(4);
    nArs = n(5);

    Ec_Td = (v/Gap_length)/Ngas*1e21;

    rates = k_Bolsig(Ec_Td);
    Te = Te_Bolsig(Ec_Td);
    mu_e_reduced = mue_N(Ec_Td);

    k1 = rates(1);
    k2 = rates(2);
    k3 = rates(3);
    k4 = rates(4);
    k5 = 8.5e-7*(Te*11600/300)^(-0.67)*1E-6;
    k6 = 6.06e-6/T*exp(-15130/T)*1E-6;
    k7 = 6.0e-10*1E-6;
    k8 = 8.75e-27*(Te)^(-4.5)*1E-12;
    k9 = 1.4e-32*1E-12;
    k10 = 2.25e-31*(T/300)^(-0.4)*1E-12;

    r1 = k1 .* nAr .* ne;           % 'Ar + e => Arp + e + e' 
    r2 = k2 .* nAr .* ne;           % 'Ar + e => Ars + e'    
    r3 = k3 .* nArs .* ne;          % 'Ars + e => Ar + e'  
    r4 = k4 .* nArs .* ne;          % 'Ars + e => Arp + e + e' 
    r5 = k5 * nAr2p .* ne;          % 'Ar2p + e => Ars + Ar'      
    r6 = k6 * nAr2p .* nAr;         % 'Ar2p + Ar => Arp + Ar + Ar' 
    r7 = k7 * nArs .* nArs;         % 'Ars + Ars => Ar2p + e'   
    r8 = k8 * nArp .* ne .* ne;     % 'Arp + e + e => Ar + e'    
    r9 = k9 * nArs .* nAr .* nAr;   % 'Ars + Ar + Ar  => Ar + Ar + Ar'     
    r10 = k10 * nArp .* nAr .* nAr; % 'Arp + Ar + Ar  => Ar2p + Ar'
    reaction_rates = [r1, r2, r3, r4, r5, r6, r7, r8, r9, r10];

    omega = reaction_rates*stoichiometric_matrix;

    I = 1.6022e-19 * ne * mu_e_reduced * Ec_Td * 1e-21 * Area;

    dydt = [omega(:); Vext - R*I - v];
end