function [ Heating ] = Heating( T_evap )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

water = importPhase('liquidVapor.xml','water');
water_vap = IdealGasMix('gri30.xml');

global open_volume porosity 
global m_solid m_metal c_metal c_solid Qst
global P_evap P_cond T_amb
global j

setState_Tsat(water,[T_evap 1]);  %vapor

set(water,'T',T_amb,'P',P_evap);
ug1 = intEnergy_mass(water);
hg1 = enthalpy_mass(water);
rhog1 = density(water);
x_gas = exergy_mass(water);
f_gas = flowExergy_mass(water);


setState_Tsat(water,[T_amb 0]); %liquid
rhoL = density(water);
rho_a1 = rhoL;

q_max = Adsorbate_Con_Ratio(T_amb,P_evap);

R = 8.314; %J/mol*K
MolarMass = 18.016;
V = open_volume;
n_init = (P_evap*V)/(R*T_amb);

m_a = q_max * m_solid;
m_gas1 = open_volume * rhog1 * porosity;
m_total = m_a + m_gas1;
n_total = (m_total/MolarMass)*1000;

P1_water = P_evap;
T1_water = T_amb;
P1_bed = P_evap;
T1_bed = T_amb;

dP = (P_cond - P_evap)/10;
i = 1;
T_prev = T1_water;
ug_prev = ug1;
ha = hg1 - Qst;  %enthalpy of adsorbate phase
ua_prev = ha - P_evap / rho_a1;
Q12 = 0;
n = n_init;

for P = P_evap:dP:P_cond
    T = T_isosteric(q_max,P);
    n = (P*V)/(R*T);
    ratio = P/T;
    
    set(water,'P',P,'T',T);
    hg = enthalpy_mass(water);
    rhog = density(water);
    ug = intEnergy_mass(water);
    sg = entropy_mass(water);
    
    setState_Tsat(water,[T 0]);  %liquid
    rhoL = density(water);
    ua = intEnergy_mass(water);
    sa = entropy_mass(water);
    
    %---------------------------------------------
    dT = T - T_prev;
    dQ_metal = m_metal * c_metal * dT;
    dQ_solid = m_solid * c_solid * dT;
    %---------------------------------------------
    dtheta = 0;
    dW = - open_volume * P * dtheta;
    
    m_a = q_max * m_solid;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P / rho_a;
    du_a = u_a - ua_prev;
    
    dQ_adsorbate = m_a * du_a - dW;
    %----------------------------------------------
    m_g = m_total - m_a;
    du_g = ug - ug_prev;
    dQ_gas = -dW + m_g * du_g;
    %-----------------------------------------------
    dQ12 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
    %-----------------------------------------------
    Q12 = Q12 + dQ12;
    
    %Exergy
    Ads.u = u_a;
    Ads.s = sa;
    Ads.v = 1/rho_a;
    Ads.m = m_a;
    Gas.u = ug;
    Gas.s = sg;
    Gas.v = 1/rhog;
    Gas.m = m_g;
    %X = exergy(T,P,Ads,Gas,dead);
    
    
    %graphing
    Heating.T_water(i) = T;
    Heating.P_water(i) = P;
    Heating.T_bed(i) = T;
    Heating.P_bed(i) = P;
    Heating.q(i) = q_max;
    Heating.m_gas(i) = m_g;
    Heating.m_ads(i) = m_a;
    Heating.Q = Q12;
    Heating.h(i) = hg;
    Heating.m_ref(j) = 0;
    Heating.P_ref(j) = P;
    
    ua_prev = u_a;
    T_prev = T;
    i = i + 1;
    j = j + 1;
end

n_final = (P*V)/(R*T);

set(water,'T',Heating.T_water(end),'P',Heating.P_water(end));

% xin = x1;
% xout = exergy_mass(water);
% fin = f1;
% fout = flowExergy_mass(water);
% 
% Xin = mdot*xin + Qg*(1+(To/Tg))
% Xout = mdot*xout

end

