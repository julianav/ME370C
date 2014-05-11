function [ T ] = T_isosteric( q,P )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Adsorption Model parameters
K = 5.5*10^-12; %Pa
Qst = 2.37 *10^6; %J/kg
R = 8314; %J/kmol*K
MolarMass = 18.016;
Rwater = R/MolarMass; %J/kg*K

T = Qst/(Rwater*log(q/(K*P)));

end

