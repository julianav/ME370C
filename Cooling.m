function [ Cooling ] = Cooling( T_max,m_init)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

water = importPhase('liquidVapor.xml','water');

global open_volume
global m_solid m_metal c_metal c_solid Qst
global P_evap P_cond To Po porosity
global j Rwater R

% set(water_vap,'P',P_evap,'T',T_amb,'X','H2O:1');
set(water,'T',T_max,'P',P_cond);
ug3 = intEnergy_mass(water);
hg3 = enthalpy_mass(water);
rhog3 = density(water);
sg3 = entropy_mass(water);

xgas3 = exergy_mass(water);
Xsolid3 = m_solid*(c_solid*(T_max-To) - To*c_solid*log(T_max/To));
Xmetal3 = m_metal*(c_metal*(T_max-To) - To*c_metal*log(T_max/To));

% setState_Psat(water,[P_evap 0]);  %liquid
setState_Tsat(water,[T_max 0]);
rhoL = density(water);
rho_a3 = rhoL;
sa3 = entropy_mass(water);

q_min = Adsorbate_Con_Ratio(T_max,P_cond);

m_a = q_min * m_solid;
m_gas = open_volume * rhog3 * porosity;
m_total = m_a + m_gas;

ha = hg3 - Qst;  %enthalpy of adsorbate phase
rho_a = rhoL;
ua_prev = ha - P_cond / rho_a;

sa3 = sg3 - (Qst/T_max-Rwater);

ref = importPhase('liquidVapor.xml','water');
set(ref,'T',To,'P',Po);
uo = intEnergy_mass(ref);
vo = 1/density(ref);
so = entropy_mass(ref);
mo = open_volume * 1/vo * porosity;

mo = m_a;
Xgas3 = m_gas*xgas3;
Xads3 = (m_a*ua_prev-mo*uo)+Po*(m_a/rho_a3*mo*vo)-To*(m_a*sa3-mo*so);
X3 = Xgas3 + Xads3 + Xsolid3 + Xmetal3;

dP = -(P_cond - P_evap)/10;
T_prev = T_max;
ug_prev = ug3;
i = 1;
Q34 = 0;
XofQ = 0;
for P = P_cond:dP:P_evap
    T = T_isosteric(q_min,P);
    
%     set(water_vap,'P',P,'T',T,'X','H2O:1');
    set(water,'P',P,'T',T);
    hg = enthalpy_mass(water);
    rhog = density(water);
    ug = intEnergy_mass(water);
    sg = entropy_mass(water);
    
    setState_Tsat(water,[T 0]);
    rhoL = density(water);
    sa = entropy_mass(water);
    
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
    
    dXofQ = dQ34*(1-To/T);
    XofQ = XofQ + dXofQ;
    
    %graphing
    Cooling.T_bed(i) = T;
    Cooling.P_bed(i) = P;
    Cooling.q(i) = q_min;
    Cooling.m_gas(i) = m_g;
    Cooling.m_ads(i) = m_a;
    Cooling.Q = Q34;
    Cooling.h(i) = hg;
%     Cooling.m_ref(j) = m_cond;
    Cooling.P_ref(j) = P;

    T_prev = T;
    ua_prev = u_a;
    ug_prev = ug;
    i = i + 1;
    j = j + 1;
end

sa = sg - (Qst/T-Rwater);

DeltaT = T-T_max;
set(water,'T',T,'P',P);
Xgas4 = m_gas*exergy_mass(water);
Xsolid4 = m_solid*(c_solid*(T-To)- To*c_solid*log(T/To));
Xmetal4 = m_metal*(c_metal*(T-To)- To*c_metal*log(T/To));
Xads4 = (m_a*u_a-mo*uo)+Po*(m_a/rho_a-mo*vo)-To*(m_a*sa-mo*so);

X4 = Xgas4 + Xads4 + Xsolid4 + Xmetal4;

DXgas = Xgas4-Xgas3;
DXsolid = Xsolid4-Xsolid3;
DXmetal = Xmetal4-Xmetal3;
DXads = Xads4-Xads3;

% fprintf('\nCooling\n')
DX = X4 - X3;
% Cooling.Xloss = Q34*(1-To/T)-DX;
Cooling.Xloss = XofQ-DX;

end

