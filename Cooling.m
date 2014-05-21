function [ Cooling ] = Cooling( T_max, m_init,m_cond )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

water = importPhase('liquidVapor.xml','water');
water_vap = IdealGasMix('gri30.xml');

global open_volume
global m_solid m_metal c_metal c_solid Qst
global P_evap P_cond
global j

% set(water_vap,'P',P_evap,'T',T_amb,'X','H2O:1');
set(water,'T',T_max,'P',P_cond);
ug3 = intEnergy_mass(water);
hg3 = enthalpy_mass(water);
rhog3 = density(water);
x3 = exergy_mass(water);

setState_Psat(water,[P_evap 0]);  %liquid
rhoL = density(water);
rho_a3 = rhoL;

q_min = Adsorbate_Con_Ratio(T_max,P_cond);

dP = -(P_cond - P_evap)/10;
T_prev = T_max;
ug_prev = ug3;
ha = hg3 - Qst;  %enthalpy of adsorbate phase
rho_a = rhoL;
ua_prev = ha - P_evap / rho_a;
i = 1;
Q34 = 0;

for P = P_cond:dP:P_evap
    T = T_isosteric(q_min,P);
    
    set(water_vap,'P',P,'T',T,'X','H2O:1');
    hg = enthalpy_mass(water_vap);
    rhog = density(water_vap);
    ug = intEnergy_mass(water_vap);
    
    setState_Tsat(water,[T 0]);
    rhoL = density(water);
    
    %---------------------------------------------------
    dT = T - T_prev;
    dQ_metal = m_metal * c_metal * dT;
    dQ_solid = m_solid * c_solid * dT;
    %---------------------------------------------------
    dtheta = 0; %Va/Vb;
    dW = - open_volume * P * dtheta;
    
    m_a = q_min * m_solid;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P / rho_a;
    du_a = u_a - ua_prev;
    
%     dQ_adsorbate = -volume * P * dtheta - m_a * du_a;
    dQ_adsorbate = m_a * du_a - dW;
    %---------------------------------------------------
    m_g = m_init;
    du_g = ug - ug_prev;
    
    dQ_gas = m_g * du_g - dW;
    %---------------------------------------------------
    dQ34 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
    %---------------------------------------------------
    Q34 = Q34 + dQ34;
    
    %graphing
    Cooling.T_bed(i) = T;
    Cooling.P_bed(i) = P;
    Cooling.q(i) = q_min;
    Cooling.m_gas(i) = m_g;
    Cooling.m_ads(i) = m_a;
    Cooling.Q = Q34;
    Cooling.h(i) = hg;
    Cooling.m_ref(j) = m_cond;
    Cooling.P_ref(j) = P;

    T_prev = T;
    ua_prev = u_a;
    ug_prev = ug;
    i = i + 1;
    j = j + 1;
end

Cooling.x_in = x3;
Cooling.x_out = exergy_mass(water);



end

