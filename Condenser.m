function [ Cond ] = Condenser(mr,T2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global To P_cond T_cond T_max

water = importPhase('liquidVapor.xml','water');

set(water,'T',T2,'P',P_cond); %incoming state
hin = enthalpy_mass(water);
FlowXin = mr*flowExergy_mass(water);

setState_Psat(water,[P_cond,1]); %vapor
Tg = temperature(water);
hg = enthalpy_mass(water);

Qsens = mr*(hg - hin);

setState_Psat(water,[P_cond,0]); %liquid
hf = enthalpy_mass(water);
FlowXout = mr*flowExergy_mass(water);

hfg = hf - hg;
Qlat = mr*hfg;

Qcond = Qsens + Qlat;

% Qcond_hfg = mr*hfg;
% Qcond = mr * (hout - hin);

XofQ = 0;
dT = (Tg-T2)/10;
for T = T2:dT:Tg
    dXofQ = Qcond/(-dT)*(1-To/T);
    XofQ = XofQ + dXofQ;
end

DX = 0;
Xloss = FlowXin - FlowXout + Qcond*(1-To/Tg) - DX;
Xloss2 = FlowXin - FlowXout + XofQ - DX;

Cond.FlowXout = FlowXout;
Cond.FlowXin = FlowXin;
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

