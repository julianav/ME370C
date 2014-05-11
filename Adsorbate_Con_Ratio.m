function[q] = Adsorbate_Con_Ratio( T, P )
%Adsorbate Consentration Ratio
%mass adsorbate / mass solid

% Adsorption Model parameters
K = 5.5*10^-12; %Pa
Qst = 2.37 *10^6; %J/kg
R = 8314; %J/kmol*K
MolarMass = 18.016;
Rwater = R/MolarMass; %J/kg*K

q = K*exp(Qst/(Rwater*T))*P;

end

