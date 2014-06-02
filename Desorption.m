function [ Desorb ] = Desorption( T_cond,T_init,T_final,m_init)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

water = importPhase('liquidVapor.xml','water');

global open_volume porosity
global m_solid m_metal c_metal c_solid Qst
global P_cond T_low P_evap
global j To Po Rwater

setState_Tsat(water,[T_cond 0]); %liquid phase
rhoL = density(water);
% sa = entropy_mass(water);

setState_Tsat(water,[T_cond,1]);
hf_low = enthalpy_mass(water);

q_max = Adsorbate_Con_Ratio(T_low,P_evap);

m_a2 = q_max * m_solid;
m_g2 = m_init;

set(water,'P',P_cond,'T',T_init);
rhog_init = density(water);

setState_Tsat(water,[T_low 0]);
rhoa_init = density(water);

%Initial state
setState_Psat(water,[P_cond,1]);
% P_cond = pressure(water);  %condensation pressure
hg2 = enthalpy_mass(water);
ug2 = intEnergy_mass(water);
sg2 = entropy_mass(water);
sa2 = sg2 - (Qst/T_init - Rwater);
set(water,'T',T_init,'P',P_cond);
Xgas2 = m_g2*exergy_mass(water);
Xsolid2 = m_solid*(c_solid*(T_init-To) - To*c_solid*log(T_init/To));
Xmetal2 = m_metal*(c_metal*(T_init-To) - To*c_metal*log(T_init/To));

ha = hg2 - Qst;  %enthalpy of adsorbate phase
rho_a = rhoL;
ua_prev = ha - P_cond / rho_a;

ref = importPhase('liquidVapor.xml','water');
set(ref,'T',To,'P',Po);
uo = intEnergy_mass(ref);
vo = 1/density(ref);
so = entropy_mass(ref);
mo = open_volume * 1/vo * porosity;
mo = m_a2;
Xads2 = (m_a2*ua_prev-mo*uo)+Po*(m_a2/rho_a-mo*vo)-To*(m_a2*sa2-mo*so);
X2 = Xgas2 + Xads2 + Xsolid2 + Xmetal2;

q_prev = q_max;
ug_prev = ug2;
ma_prev = m_a2;
mg_prev = m_g2;
Qcond = 0;
Q23 = 0;
Qgas = 0;
Qads = 0;
XofQ = 0;
m_des = 0;
m_leaving_prev = 0;
m_cond = 0;
FlowXin = 0;
FlowXout = 0;
i = 1;
T_prev = T_init;
theta_prev = (m_a2/rhoL)/open_volume;
dT = (T_final - T_init)/10;
% T_final - T_init
for T = T_init:dT:T_final
    P = P_cond;
    set(water,'T',T,'P',P);
%     setState_Psat(water,[P_cond 1]);
    hg = enthalpy_mass(water);
    ug = intEnergy_mass(water);
    rhog = density(water);
    sg = entropy_mass(water);   
    
    set(water,'T',T,'P',P);
    flow_x = flowExergy_mass(water);
        
    mg_capacity = open_volume * rhog * porosity;

    setState_Tsat(water,[T 0]); %liquid
    rhoL = density(water);
    sa = entropy_mass(water);
    %-----------------------------------------
    dT = T - T_prev;
    dQ_metal = m_metal * c_metal * dT;
    dQ_solid = m_solid * c_solid * dT;
    %-----------------------------------------
    q = Adsorbate_Con_Ratio(T,P);
    Dq = q - q_prev;
    
    m_a = q * m_solid;
%     Dm_a = m_a - ma_prev;
Dm_a = m_solid*Dq;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P_cond / rho_a;
    du_a = u_a - ua_prev;
    dU_a = m_a * du_a + u_a * Dm_a;
    
    Dm_des = -Dm_a;
     
    theta = (m_a/rhoL)/open_volume;  
    dtheta = theta - theta_prev;
    dW_a = open_volume * P_cond * dtheta;%Dq;

    dQdes = Dm_des*Qst;
    
%     dQ_adsorbate = dU_a - dW_a - Dm_a * hg;  %changed mass of des to mass of ads
       dQ_adsorbate = dU_a + dW_a + Dm_des * hg;  %changed mass of des to mass of ads
    Qads = Qads + dQ_adsorbate;
       %---------------------------------------------
    m_g = open_volume*rhog_init*(porosity+(m_solid*q)/(open_volume*rho_a));
    
    Dm_g = m_g - mg_prev;
    
Dm_cond = -(Dm_a + Dm_g);
dqcond = hg * Dm_cond;
    

%     mg_init = mg_prev + Dm_g; % mass of gas initially
%     
%     if mg_init > mg_capacity;   % excess mass leaves the adsorber
%         m_leaving = mg_init - mg_capacity;
%     else
%         m_leaving = 0;
%     end
%     
%     m_g = mg_init - m_leaving;
%     Dm_g = m_g - mg_prev;
    
    du_g = ug - ug_prev;
    dU_g = m_g * du_g + ug * Dm_g;  %m_g = Dm_a - m_leaving
    
    dW_g = -dW_a;
    
%     Dm_cond = m_leaving;
    m_cond = m_cond + Dm_cond;
    dQ_cond = Dm_cond * (hg - hf_low);
    
    dQ_gas = dU_g - dW_g - Dm_g * hg;
    
    Qgas = Qgas + dQ_gas;
    %---------------------------------------------
    dQ23 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
    %---------------------------------------------

    %---------------------------------------------
    Q23 = Q23 + dQ23;   
    Qcond = Qcond + dQ_cond;
    
    m_des = m_des + Dm_des;
    
    m_leaving = Dm_cond;
    dFlowXin = flow_x * Dm_g;
    dFlowXout = flow_x * m_leaving;
    
    FlowXin = FlowXin + dFlowXin;
    FlowXout = FlowXout + dFlowXout;
    
    dXofQ = dQ23*(1-To/T);
    XofQ = XofQ + dXofQ;  
    
    
    %graphing
%     Desorb.T_water(i) = T;  %?????
    Desorb.P_water(i) = P_cond;
    Desorb.T_bed(i) = T;
    Desorb.P_bed(i) = P_cond;
    Desorb.q(i) = q;
    Desorb.m(i) = m_des;
    Desorb.m_gas(i) = m_g;
    Desorb.m_ads(i) = m_a;
    Desorb.m_cond = m_cond;  
    Desorb.Q_gas(i) = Qgas;
    Desorb.h(i) = hg;
    Desorb.m_ref = m_cond;
    Desorb.P_ref(j) = P_cond;
    Desorb.Q_cond = Qcond;
    Desorb.Q_des = Q23;
    
    T_prev = T;
    ua_prev = u_a;
    ug_prev = ug;
    q_prev = q;
    ma_prev = m_a;
    mg_prev = m_g;
    m_leaving_prev = m_leaving;
    theta_prev = theta;
    
    i = i + 1;
    j = j + 1;
end

% mgss = max(Desorb.m_gas)

sa = sg - (Qst/T - Rwater);
mo = m_a;
set(water,'T',T_final,'P',P_cond);
Xgas3 = m_g*exergy_mass(water);
Xsolid3 = m_solid*(c_solid*(T_final-To) - To*c_solid*log(T_final/To));
Xmetal3 = m_metal*(c_metal*(T_final-To) - To*c_metal*log(T_final/To));
Xads3 = (m_a*u_a-mo*uo)+Po*(m_a/rho_a-mo*vo)-To*(m_a*sa-mo*so);

X3 = Xgas3 + Xads3 + Xsolid3 + Xmetal3;

% fprintf('\nDesorb\n')
DX = X3-X2;

%Heating stream
% set(water,'T',Thotstream,'P',Po);
% heating_hin = enthalpy_mass(water);
% 
% heating_hout = heating_hin - Q23/mdot_heating;
% set(water,'P',Po,'H',heating_hout);
% heatingToutFinal = temperature(water)-273;


% Desorb.Xloss = Q23*(1-To/T_final) - FlowXout - DX;
Desorb.Xloss = XofQ - FlowXout - DX;
Desorb.DX = DX;
Desorb.FlowXout = FlowXout;
end

% set(water,'T',T_cond,'P',P_cond);
% xin = x_init;
% xout = exergy_mass(water);
% 
% set(water,'T',T_cond,'P',P_cond);
% fin = f_init;
% fout = flowExergy_mass(water);
% 
% mdot = max(Desorb.m_ref);
% Qg = Qgas;
% Xin = mdot*fin + Qg*(1-(To/T_init));
% Xout = mdot*fout;
% 
% Desorb.xin = fin;
% Desorb.xout = fout;
% Desorb.Xloss = Xin - Xout;

