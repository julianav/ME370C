function [ Exp ] = Expansion( m )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global P_evap T_evap T_cond P_cond

water = importPhase('liquidVapor.xml','water');

%incoming state

% T3_water = T_cond;
% P3_water = P_cond;
% set(water,'T',T3_water,'P',P3_water);
% h3 = enthalpy_mass(water);
% h = h3;

setState_Psat(water,[P_cond,0]); %liquid
hin = enthalpy_mass(water);
FlowXin = m*flowExergy_mass(water);


%exiting state
hout = hin;
set(water,'P',P_evap,'H',hin);
FlowXout = m*flowExergy_mass(water);

setState_Psat(water,[P_evap,0]); %liquid
hf = enthalpy_mass(water);
setState_Psat(water,[P_evap,1]); %vapor
hg = enthalpy_mass(water);
T = temperature(water);
P = P_evap;
x = (hout - hf)/(hg - hf);
y = 1-x;

hfg = hg - hf;

Xloss = FlowXin - FlowXout;

Exp.T = T;
Exp.P = P;
Exp.y = y;
Exp.hfg = hfg;
Exp.hout = hout;
% Exp.xout = FlowXout/m;
Exp.XexpDest = Xloss;

end

