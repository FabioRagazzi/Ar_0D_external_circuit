clearvars, close all, clc %#ok<DUALC>

eps0 = 8.854e-12;
e = 1.6022e-19;
kB = 1.38e-23;

Gap_length = 4e-3;
Radius = 4e-3;
Vext = 1e3;
R = 1e5;
T = 300;
Ngas = 101325*100/760/(kB*T);
Area = pi * Radius^2;

Te_Bolsig = load('ZDP_m_c_mean_en.txt'); 
Te_Bolsig = griddedInterpolant(Te_Bolsig(:,1),Te_Bolsig(:,2)*2/3,"linear","nearest");

mue_N = load('ZDP_m_c_mue_N.txt'); 
mue_N = griddedInterpolant(mue_N(:,1),mue_N(:,2),"pchip","nearest");

k_Ar_iz = load('ZDP_m_c_k_Ar_iz.txt');
k_Ar_exc = load('ZDP_m_c_k_Ar_exc.txt'); 
k_Ars_de_exc = load('ZDP_m_c_k_Ars_de_exc.txt'); 
k_Ars_iz = load('ZDP_m_c_k_Ars_iz.txt');
k_Bolsig = griddedInterpolant(k_Ars_iz(:,1),[k_Ar_iz(:,2) k_Ar_exc(:,2) k_Ars_de_exc(:,2) k_Ars_iz(:,2)]); % m3/s rate Coefficients

% reactions
% 'Ar + e -> Arp + e + e';    
% 'Ar + e -> Ars + e';        
% 'Ars + e -> Ar + e';       
% 'Ars + e -> Arp + e + e';  
% 'Ar2p + e -> Ars + Ar';         
% 'Ar2p + Ar -> Arp + Ar + Ar';   
% 'Ars + Ars -> Ar2p + e';        
% 'Arp + e + e -> Ar + e';        
% 'Ars + Ar + Ar -> Ar + Ar + Ar';
% 'Arp + Ar + Ar -> Ar2p + Ar';   
stoichiometric_matrix = [1    -1     0     1     0;
                         0    -1     0     0     1;
                         0     1     0     0    -1;
                         1     0     0     1    -1;
                        -1     1    -1     0     1;
                         0     1    -1     1     0;
                         1     0     1     0    -2;
                        -1     1     0    -1     0;
                         0     1     0     0    -1;
                         0    -1     1    -1     0;];

%    [e,    Ar,   Ar2p, Arp, Ars]
n0 = [1e12, Ngas, 1e9, 1e12, 1e12]';
y0 = [n0(:); Vext];

Mass = eye(6);
Mass(end,end) = 0;
opts = odeset('Mass',Mass);

t_end = 1e-3;
fode = @(t,y)OdeFuncChemistry(t,y,k_Bolsig,stoichiometric_matrix,Te_Bolsig,mue_N,T,Ngas,Gap_length,Vext,R,Area);
[tout,yout] = ode15s(fode,linspace(0,t_end,2), y0(:), opts);

% plot results
% --------------------------------------------------------------------------------------------------
nout = yout(:,1:5) * 1e-6;
vout = yout(:,end);
ETdout = (vout/Gap_length)/Ngas*1e21;

fig = figure;
ax = axes(fig);
colors = ax.ColorOrder;
yyaxis left
loglog(tout,nout(:,1),"LineWidth",2,"Color",colors(1,:),"LineStyle","-","DisplayName","e")
hold on
loglog(tout,nout(:,3),"LineWidth",2,"Color",colors(4,:),"LineStyle","-","DisplayName","Ar^{2+}")
loglog(tout,nout(:,4),"LineWidth",2,"Color",colors(3,:),"LineStyle","-","DisplayName","Ar^{+}")
loglog(tout,nout(:,5),"LineWidth",2,"Color",colors(2,:),"LineStyle","-","DisplayName","Ar^{*}")
ylabel("Number density (cm^{-3})")
yyaxis right
semilogx(tout,ETdout,"LineWidth",2,"Color","green","LineStyle","--","Color",colors(5,:),"DisplayName","E")
ylabel("Reduced electric field (Td)")
legend("Location","bestoutside")
grid on
ax.FontSize = 18;
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = colors(5,:);
xlim([1e-9,t_end])

