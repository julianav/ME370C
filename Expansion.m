function [ Exp ] = Expansion( h,m,xin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global P_evap

water = importPhase('liquidVapor.xml','water');

setState_Psat(water,[P_evap,0]); %liquid
hf = enthalpy_mass(water);
setState_Psat(water,[P_evap,1]); %vapor
hg = enthalpy_mass(water);
T = temperature(water);
P = P_evap;
x = (h - hf)/(hg - hf);
y = 1-x;

hfg = hg - hf;

set(water,'P',P_evap,'H',h);
xout = exergy_mass(water);

Xin = m*xin;
Xout = m*xout;
Xloss = Xin - Xout;


Exp.T = T;
Exp.P = P;
Exp.y = y;
Exp.hfg = hfg;
Exp.xout = xout;
Exp.Xloss = Xloss;

end

