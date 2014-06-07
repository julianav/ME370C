function [ Heating ] = Heating( hotTin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

water = importPhase('liquidVapor.xml','water');
% water_vap = IdealGasMix('gri30.xml');

global open_volume porosity 
global m_solid m_metal c_metal c_solid Qst
global P_evap P_cond T_low To Po
global j Rwater R 

% setState_Tsat(water,[T_evap 1]); 

set(water,'T',T_low,'P',P_evap);
ug1 = intEnergy_mass(water);
hg1 = enthalpy_mass(water);
rhog1 = density(water);
sg1 = entropy_mass(water);
sa1 = sg1 - (Qst/T_low - Rwater);

xgas1 = exergy_mass(water);
Xsolid1 = m_solid*(c_solid*(T_low-To) - To*c_solid*log(T_low/To));
Xmetal1 = m_metal*(c_metal*(T_low-To) - To*c_metal*log(T_low/To));

setState_Tsat(water,[T_low 0]); %liquid
% Ptest = pressure(water);
rhoL = density(water);
rho_a1 = rhoL;
% sa1 = entropy_mass(water);
% uatest = intEnergy_mass(water);
% xadstest = exergy_mass(water);

q_max = Adsorbate_Con_Ratio(T_low,P_evap);

R1 = R/1000; %J/mol*K
MolarMass = 18.016;
V = open_volume;
n_init = (P_evap*V)/(R1*T_low);

m_a = q_max * m_solid;
m_gas = open_volume * rhog1 * porosity;
m_total = m_a + m_gas;
n_total = (m_total/MolarMass)*1000;

% P1_water = P_evap;
% T1_water = T_low;
% P1_bed = P_evap;
% T1_bed = T_low;

ug_prev = ug1;
ha = hg1 - Qst;  %enthalpy of adsorbate phase
ua_prev = ha - P_evap / rho_a1;

ref = importPhase('liquidVapor.xml','water');
set(ref,'T',To,'P',Po);
uo = intEnergy_mass(ref);
vo = 1/density(ref);
so = entropy_mass(ref);
mo = open_volume * 1/vo * porosity;
mo = m_a;
Xgas1 = m_gas*xgas1;
Xads1 = (m_a*ua_prev-mo*uo)+Po*(m_a/rho_a1-mo*vo)-To*(m_a*sa1-mo*so);
X1 = Xgas1 + Xads1 + Xsolid1 + Xmetal1;

T_prev = T_low;
vg_prev = 1/rhog1;
theta_prev = (m_a/rhoL)/open_volume;
dP = (P_cond - P_evap)/10;
Q12 = 0;
XofQ = 0;
n = n_init;
i = 1;
for P = P_evap:dP:P_cond
    T = T_isosteric(q_max,P);
    n = (P*V)/(R*T);
    ratio = P/T;
    
    set(water,'P',P,'T',T);
    setState_Tsat(water,[T,1]);
    hg = enthalpy_mass(water);
    vg = 1/density(water);
    ug = intEnergy_mass(water);
    sg = entropy_mass(water);
%     xgas = exergy_mass(water);
    
    setState_Tsat(water,[T 0]);  %liquid
    rhoL = density(water);
%     ua = intEnergy_mass(water);
    sa = entropy_mass(water);
    
    %---------------------------------------------
    dT = T - T_prev;
    dQ_metal = m_metal * c_metal * dT;
    dQ_solid = m_solid * c_solid * dT;
    %---------------------------------------------
    theta = (m_a/rhoL)/open_volume;
    dtheta = theta - theta_prev;
    dW = - open_volume * P * dtheta;
    
    m_a = q_max * m_solid;
    
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P / rho_a;
    du_a = u_a - ua_prev;
    
%     Du = (u_a-uo);
%     Dv = (1/rho_a-vo);
%     Ds = (sa-so);
%     xads2 = (u_a-uo)+Po*(1/rho_a-vo)-To*(sa-so);
%    
    dQ_adsorbate = m_a * du_a - dW;
    %----------------------------------------------
%     m_g = m_total - m_a;
    du_g = ug - ug_prev;
    dQ_gas = -dW + m_gas * du_g;
    %-----------------------------------------------
    dQ12 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
    %-----------------------------------------------
    Q12 = Q12 + dQ12;
    
    dXofQ = dQ12*(1-To/T);
    XofQ = XofQ + dXofQ;  
    
    %graphing
    Heating.T_water(i) = T;
    Heating.P_water(i) = P;
    Heating.T_bed(i) = T;
    Heating.P_bed(i) = P;
    Heating.q(i) = q_max;
    Heating.m_gas(i) = m_gas;
    Heating.m_ads(i) = m_a;
    Heating.Q = Q12;
    Heating.h(i) = hg;
    Heating.m_ref(j) = 0;
    Heating.P_ref(j) = P;
    
    ua_prev = u_a;
    T_prev = T;
    theta_prev = theta;
    i = i + 1;
    j = j + 1;
end

n_final = (P*V)/(R1*T);

sa = sg - (Qst/T-Rwater);

T_final = Heating.T_water(end);
set(water,'T',T,'P',P);
Xgas2 = m_gas*exergy_mass(water);
Xsolid2 = m_solid*(c_solid*(T-To)- To*c_solid*log(T/To));
Xmetal2 = m_metal*(c_metal*(T-To)- To*c_metal*log(T/To));
Xads2 = (m_a*u_a-mo*uo)+Po*(m_a/rho_a-mo*vo)-To*(m_a*sa-mo*so);

X2 = Xgas2 + Xads2 + Xsolid2 + Xmetal2;

DXgas = Xgas2-Xgas1;
DXsolid = Xsolid2-Xsolid1;
DXmetal = Xmetal2-Xmetal1;
DXads = Xads2-Xads1;


%Heating stream
% set(water,'T',hotTin,'P',Po);
% heating_hin = enthalpy_mass(water);
% Qheat = -Q12;
% heating_hout = heating_hin - Q12/mdot_heating;
% set(water,'P',Po,'H',heating_hout);
% heatingTout = temperature(water);
% heatout = heatingTout-273;




% fprintf('\nHeating\n')
DX = X2 - X1;

% Heating.Xloss = Q12*(1-To/T) - DX;
Heating.Xloss = XofQ - DX;
Heating.DX = DX;
% Heating.heatingTout = heatingTout;

end

% xin = x1;
% xout = exergy_mass(water);
% fin = f1;
% fout = flowExergy_mass(water);
% 
% Xin = mdot*xin + Qg*(1+(To/Tg))
% Xout = mdot*xout
