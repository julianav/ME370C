function[q] = Adsorbate_Con_Ratio( T, P )
%Adsorbate Consentration Ratio
%mass adsorbate / mass solid

global K Qst

% Adsorption Model parameters
% Qst = 2.51 *10^6;
R1 = 8314; %J/kmol*K
MolarMass = 18.016;
Rwater = R1/MolarMass; %J/kg*K

q_old = K*exp(Qst/(Rwater*T))*P;
q = q_old;

% Ko = 7.3e-10; %kPa-1
% Q = 2693; %kJ/kg 
% R = 0.461;
% q_eq = 0.4;
% t = 6;
% y = exp(Q/(R*T));
% 
% q = (Ko*y*P/1e3)/(1+((Ko/q_eq)*y*P/1e3)^t)^(1/t);



end

