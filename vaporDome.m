function [h,s,u,T,P,v] = vaporDome(species);
%    This function creates a vapor dome of the input species
%    specified by the user. The function is called as follows:
%
%    [domeProps] = VaporDome(species)
%
%    where "species" is the species name of the pure substance
%
%    To plot Vapor Dome on P-V diagram simply type:
%    plot(v,P)
%
%    Paul Mobley
%    Winter 2010
%
%    WARNING: This specific function is currently coded for a CO2
%    vapor dome.  To use for other substances you may have to change line
%    where T_min is set.

T_crit = critTemperature(species);
T_min = minTemp(species)+17;
T_vector = linspace(T_min, T_crit-0.1);
i=1;

%Liquid Line
for T_liq = T_vector
     Liquid_dome = setState_Tsat(species,[T_liq 0]);
     h_liq(i,1) = enthalpy_mass(Liquid_dome);
     s_liq(i,1) = entropy_mass(Liquid_dome);
     u_liq(i,1) = intEnergy_mass(Liquid_dome);
     v_liq(i,1) = 1/density(Liquid_dome);
     P_liq(i,1) = pressure(Liquid_dome);
     i=i+1;
end
i=1;

%Vapor Line
for T_vap = T_vector
     Vapor_dome = setState_Tsat(species,[T_vap 1]);
     h_vap(i,1) = enthalpy_mass(Vapor_dome);
     s_vap(i,1) = entropy_mass(Vapor_dome);
     u_vap(i,1) = intEnergy_mass(Vapor_dome);
     v_vap(i,1) = 1/density(Vapor_dome);     
     P_vap(i,1) = pressure(Vapor_dome);
     i=i+1;
end

% Vapor Dome values that can be plotted against each other
T_vector = rot90(T_vector,-1);
h = [h_liq; flipud(h_vap)];
s = [s_liq; flipud(s_vap)];
u = [u_liq; flipud(u_vap)];
T = [T_vector; flipud(T_vector)];
P = [P_liq; flipud(P_vap)];
v = [v_liq; flipud(v_vap)];
