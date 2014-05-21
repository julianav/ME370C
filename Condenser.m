function [ Cond ] = Condenser(mr,mc,Tc,xin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global Po P_cond

water = importPhase('liquidVapor.xml','water');

setState_Psat(water,[P_cond,0]); %liquid
hf = enthalpy_mass(water);
xout = exergy_mass(water);
setState_Psat(water,[P_cond,1]); %vapor
hg = enthalpy_mass(water);
hfg = hg - hf;

Qcond = mr * hfg;

%Condenser
cond_water = importPhase('liquidVapor.xml','water');
set(cond_water,'T',Tc,'P',Po);
hin = enthalpy_mass(cond_water);
xcin = exergy_mass(cond_water);
hout = hin + (Qcond/1e3)/mc;
set(cond_water,'P',Po,'H',hout);
Tout = temperature(cond_water);
xcout = exergy_mass(cond_water);

Xin = mr*xin + mc*xcin;
Xout = mr*xout + mc*xcout;
Xloss = Xin - Xout;

Cond.xout = xout;
Cond.Qcond = Qcond;
Cond.Xloss = Xloss;

end

