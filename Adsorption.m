function [ Adsorb ] = Adsorption( T_init,T_final,m_init,m_cond )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

water = importPhase('liquidVapor.xml','water');
water_vap = IdealGasMix('gri30.xml');

global open_volume porosity
global m_solid m_metal c_metal c_solid Qst
global P_cond P_evap T_max
global j Po To

T4_bed = T_init;
set(water,'T',T_init,'P',P_evap);
% set(water_vap,'P',P_evap,'T',T_init,'X','H2O:1');
ug4 = intEnergy_mass(water);
hg4 = enthalpy_mass(water);
rhog4 = density(water);

xgas4 = exergy_mass(water);
Xsolid4 = m_solid*(c_solid*(T_init-To) - To*c_solid*log(T_init/To));
Xmetal4 = m_metal*(c_metal*(T_init-To) - To*c_metal*log(T_init/To));

setState_Tsat(water,[T4_bed 0]); %liquid phase
rhoL = density(water);
sa = entropy_mass(water);

q_min = Adsorbate_Con_Ratio(T_max,P_cond);

m_a4 = q_min * m_solid;
m_g4 = m_init;

ha = hg4 - Qst;  %enthalpy of adsorbate phase
rho_a = rhoL;
ua_prev = ha - P_cond / rho_a;

ref = importPhase('liquidVapor.xml','water');
set(ref,'T',To,'P',Po);
uo = intEnergy_mass(ref);
vo = 1/density(ref);
so = entropy_mass(ref);
mo = open_volume * 1/vo * porosity;

mo=m_a4;
Xads4 = (m_a4*ua_prev-mo*uo)+Po*(m_a4/rho_a-mo*vo)-To*(m_a4*sa-mo*so);
X4 = m_g4*xgas4 + Xads4 + Xsolid4 + Xmetal4;

q_prev = q_min;
ug_prev = ug4;
ma_prev = m_a4;
mg_prev = m_g4;
T_prev = T4_bed;
Q41 = 0;
m_ads = 0;
m_evap = 0;
Q_evap = 0;
Qsens = 0;
Qgas = 0;
FlowXin = 0;
FlowXout = 0;
i = 1;
dT = -(T_init - T_final)/10;

for T = T_init:dT:T_final
    q = Adsorbate_Con_Ratio(T,P_evap);
    
    set(water,'T',T,'P',P_evap);
    hg = enthalpy_mass(water);
    rhog = density(water);
    ug = intEnergy_mass(water);
    flow_x = flowExergy_mass(water);

    mg_capacity = open_volume * rhog * porosity;
    
    setState_Psat(water,[P_evap,1]);
    h_ev = enthalpy_mass(water);
    
    setState_Tsat(water,[T 0]);
    rhoL = density(water);
    
    %--------------------------------------------------
%     dT = T - T_prev;
    dQ_metal = m_metal * c_metal * dT;
    dQ_solid = m_solid * c_solid * dT;
    %--------------------------------------------------
    Dq = q - q_prev;
    m_a = q * m_solid;
    Dm_a = m_a - ma_prev;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P_evap / rho_a;
    du_a = u_a - ua_prev;
    dU_a = m_a * du_a + u_a * Dm_a;

    Dm_ads = Dm_a;
    
    dW_a = - open_volume * P_evap * Dq;
    
    dQ_adsorbate = -h_ev * Dm_ads + dU_a - dW_a;
    %-------------------------------------------------
%     m_g_paper = (volume * density_amb) * (porosity - (q_max * m_solid)/(volume * rho_a1));

   
    m_g = mg_capacity;
    Dm_g = m_g - mg_prev;
    Dm_entering = Dm_ads - Dm_g;
 
    du_g = ug - ug_prev;
    
    dU_g = m_g * du_g + ug * Dm_g;
    
    dQ_gas = dU_g - dW_a - Dm_g * h_ev;
    Qgas = Qgas + dQ_gas;
    %--------------------------------------------------
    dQ41 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
    %--------------------------------------------------
    Dm_evap = Dm_entering;
    m_evap = m_evap + Dm_evap;
    
    dQ_evap = Dm_evap * (h_ev);
    %--------------------------------------------------
    Q41 = Q41 + dQ41;
    Q_evap = Q_evap + dQ_evap;
    
    Qsens = Qsens + dQ_metal + dQ_solid;
    
    m_ads = m_ads + Dm_ads;
    
    dFlowXin = flow_x * Dm_entering;
    dFlowXout = flow_x * Dm_ads;
    
    FlowXin = FlowXin + dFlowXin;
    FlowXout = FlowXout + dFlowXout;

    
    %graphing
%     T41_water(i) = T;  %?????
    Adsorb.P_water(i) = P_evap;
    Adsorb.T_bed(i) = T;
    Adsorb.P_bed(i) = P_evap;
    Adsorb.q(i) = q;
    Adsorb.m_ads_tot(i) = m_ads;
    Adsorb.m_gas(i) = m_g;
    Adsorb.m_ads(i) = m_a;
    Adsorb.m_evap(i) = m_evap;
    Adsorb.Q = Q41;
    Adsorb.h(i) = hg;
    Adsorb.m_ref(j) = m_cond - m_evap;
    Adsorb.P_ref(j) = P_evap;
    Adsorb.Q_sens = Qsens;
  
    T_prev = T;
    q_prev = q;
    ua_prev = u_a;
    ug_prev = ug;
    ma_prev = m_a;
    mg_prev = m_g;
    
    i = i + 1;
    j = j + 1;
end

set(water,'T',T,'P',P_evap);
mo = m_a;
Xgas1 = m_g*exergy_mass(water);
Xsolid1 = m_solid*(c_solid*(T_final-To) - To*c_solid*log(T_final/To));
Xmetal1 = m_metal*(c_metal*(T_final-To) - To*c_metal*log(T_final/To));
Xads1 = (m_a*u_a-mo*uo)+Po*(m_a/rho_a-mo*vo)-To*(m_a*sa-mo*so);

X1 = Xgas1 + Xads1 + Xsolid1 + Xmetal1;

DX = X1 - X4;
% DFlowX = FlowXin - FlowXout;

Adsorb.Xloss = FlowXin - DX;

end


% 
% set(water,'T',T,'P',P_evap);
% fin = f_init;
% fout = flowExergy_mass(water);
% 
% mdot = max(Adsorb.m_ref);
% Qg = Qgas;
% Xin = mdot*fin + Qg*(1-(To/T_init));
% Xout = mdot*fout;
% 
% Adsorb.xin = fin;
% Adsorb.xout = fout;
% Adsorb.Xloss = Xin - Xout;