function [ Cond ] = Condenser(mr,mc,T2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global To P_cond T_cond

water = importPhase('liquidVapor.xml','water');

set(water,'T',T2,'P',P_cond); %incoming state
hin = enthalpy_mass(water);
xin = flowExergy_mass(water);

setState_Psat(water,[P_cond,1]); %vapor
Tg = temperature(water);
hg = enthalpy_mass(water);

Qsens = mr*(hg - hin);

setState_Psat(water,[P_cond,0]); %liquid
hf = enthalpy_mass(water);
xout = flowExergy_mass(water);

hfg = hf - hg;
Qlat = mr*hfg;

Qcond = Qsens + Qlat;

% Qcond_hfg = mr*hfg;
% Qcond = mr * (hout - hin);

Xin = mr*xin;
Xout = mr*xout;
Xloss = Xin - Xout;% + Qcond*(1-To/T_cond);

Cond.xout = xout;
Cond.Qcond = Qcond;
% Cond.Qcond_hfg = Qcond_hfg;
Cond.Xloss = Xloss;

end

% %Condenser
% cond_water = importPhase('liquidVapor.xml','water');
% set(cond_water,'T',Tc,'P',Po);
% hin = enthalpy_mass(cond_water);
% xcin = exergy_mass(cond_water);
% hout = hin + (Qcond/1e3)/mc;
% set(cond_water,'P',Po,'H',hout);
% Tout = temperature(cond_water);
% xcout = exergy_mass(cond_water);

