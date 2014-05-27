function [ Evap ] = Evaporator(hin,mr,T1)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

global To P_evap P_cond T_evap T_cond 

water = importPhase('liquidVapor.xml','water');

%incoming state
set(water,'P',P_evap,'H',hin);
y = (1-vaporFraction(water));
xin = flowExergy_mass(water);

setState_Psat(water,[P_evap,0]);
hf = enthalpy_mass(water);
setState_Psat(water,[P_evap,1]);
hg = enthalpy_mass(water);
hfg = hf - hg;

Qlat = -mr*hfg*y;

set(water,'T',T1,'P',P_evap);
hout = enthalpy_mass(water);
xout = flowExergy_mass(water);

Qsens = mr*(hout - hg);

Qevap = Qlat + Qsens;

Xin = mr*xin ;
Xout = mr*xout;
Xloss = Xin - Xout;% + Qevap*(1-(To/T_evap));



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








%%
Evap.Qevap = Qevap;
% Evap.Qevap_hfg = Qevap_hfg;
Evap.Xloss = Xloss;


end

% q_max = Adsorbate_Con_Ratio(T_amb,P_evap);
% q_min = Adsorbate_Con_Ratio(T_max,P_cond);
% Qevap = (q_max - q_min) * hfg * m_solid ;
