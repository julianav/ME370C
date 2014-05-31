function [ Evap ] = Evaporator(hin,mr,T1)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

global To P_evap P_cond T_evap T_cond 

water = importPhase('liquidVapor.xml','water');

%incoming state
set(water,'P',P_evap,'H',hin);
y = (1-vaporFraction(water));
FlowXin = mr*flowExergy_mass(water);
h1 = enthalpy_mass(water);
s1 = entropy_mass(water);

% setState_Psat(water,[P_evap,0]);
% hf = enthalpy_mass(water);
setState_Psat(water,[P_evap,1]);
hg = enthalpy_mass(water);
% hfg = hf - hg;

% Qlat = -mr*hfg*y;
Qlat = mr*(hg-hin);

set(water,'T',T1,'P',P_evap);
% setState_Psat(water,[P_evap 1]);
hout = enthalpy_mass(water);
FlowXout = mr*flowExergy_mass(water);
h2 = enthalpy_mass(water);
s2 = entropy_mass(water);

Qsens = mr*(hout - hg);

Qevap = Qlat + Qsens;

DX = 0;

% dFlowX = mr*((h1-h2)-To*(s1-s2))
% dFlowX2 = FlowXin - FlowXout

Xloss = FlowXin + Qevap*(1-(To/T_evap)) - FlowXout - DX;

Evap.Qevap = Qevap;
Evap.Xloss = Xloss;
Evap.Flowxout = FlowXout/mr;
Evap.FlowXout = FlowXout;
Evap.FlowXin = FlowXin;

%%

% set(water,'T',T_evap,'P',P_evap);
% xout = flowExergy_mass(water);
% 
% setState_Tsat(water,[T_evap 1]);
% hge = enthalpy_mass(water);
% sge = entropy_mass(water);
% xge = flowExergy_mass(water);
% 
% setState_Tsat(water,[T_cond 0]);
% hf_init = enthalpy_mass(water);
% sf_init = entropy_mass(water);
% setState_Tsat(water,[T_evap 0]);
% hf_flash = enthalpy_mass(water);
% sf_flash = entropy_mass(water);
% 
% mr_init = mr; 
% mr_flash = mr*y;
% Xexp = (mr_init - mr_flash)*hfg*(To/T_evap - 1);
% Xloss_exp = -(mr_init-mr_flash)*(hge - To*sge)...
%     - (mr_flash*(hf_flash-To*sf_flash)-mr_init*(hf_init-To*sf_init));
% 
% Xevap = mr_flash*hfg*(To/T_evap - 1);










end

% q_max = Adsorbate_Con_Ratio(T_amb,P_evap);
% q_min = Adsorbate_Con_Ratio(T_max,P_cond);
% Qevap = (q_max - q_min) * hfg * m_solid ;
