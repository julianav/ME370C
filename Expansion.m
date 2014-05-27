function [ Exp ] = Expansion( h,m )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global P_evap T_evap T_cond P_cond

water = importPhase('liquidVapor.xml','water');

%incoming state
setState_Psat(water,[P_cond,0]); %liquid
hin = enthalpy_mass(water);
xin = flowExergy_mass(water);


%exiting state
hout = hin;
set(water,'P',P_evap,'H',h);
xout = flowExergy_mass(water);

setState_Psat(water,[P_evap,0]); %liquid
hf = enthalpy_mass(water);
setState_Psat(water,[P_evap,1]); %vapor
hg = enthalpy_mass(water);
T = temperature(water);
P = P_evap;
x = (h - hf)/(hg - hf);
y = 1-x;

hfg = hg - hf;

Xin = m*xin;
Xout = m*xout;
Xloss = Xin - Xout;

Exp.T = T;
Exp.P = P;
Exp.y = y;
Exp.hfg = hfg;
Exp.hout = h;
Exp.xout = xout;
Exp.Xloss = Xloss;

end

