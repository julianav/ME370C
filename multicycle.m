function[Outputs] = multicycle(hot_T_in,chilled_T_in,cooling_T_in,hot_T_diff,chilled_T_diff,cooling_T_diff)

% clear all
% format compact
% fprintf('\n***************************************************************\n')
% fprintf('NEW RUN')
% fprintf('\n***************************************************************\n')

water = importPhase('liquidVapor.xml','water');
% water_vap = IdealGasMix('gri30.xml');

global open_volume porosity 
global m_solid m_metal c_metal c_solid
global P_evap P_cond T_low T_max T_evap T_cond
global j cycletime
global To Po
% global dead
global Qst K Rwater R

T_low = 300;
To = 298;
Po = 101325;
% set(water_vap,'T',dead.To,'P',dead.Po,'X','H2O');
% set(water,'T',dead.To,'P',dead.Po);
% dead.mu_o = chemPotentials(water);
% dead.so = entropy_mass(water);
% dead.uo = intEnergy_mass(water);
% dead.vo = 1/density(water);

cycletime = 6*60; %second %3minutes

%Properties
volume = 2.5; %3.1*.9*.9;  %m^3 chamber volume find value
adsorbent_density = 800; %800kg/m^3
porus_volume = 0.35; %ml/kg
water_content = 0; %percent
porosity = adsorbent_density*porus_volume*1000*1e-6; %fraction

%Masses
Ntubes = 672;
V_tubes = 3.4 * pi * ((16e-3/2)^2 - ((16e-3/2) - .8e-3)^2) * Ntubes;
density_Cu = 8940; %kg/m^3
Nfins = 1.2e5;
V_fins = 340e-3 * 28e-3 * .95e-3 * Nfins;
density_Al = 2700;

m_Cu = V_tubes * density_Cu;
m_Al = V_fins * density_Al;
m_metal = m_Cu + m_Al;
volume_displace = pi * (16e-3/2)^2 * 3.4 * Ntubes + V_fins;
open_volume = volume - volume_displace;
m_solid = (open_volume) * adsorbent_density; %kg    % silica_Mass = 895; %kg
c_Cu = 384.4; %J/kg*K  %metal tubes, copper
c_Al = 904; %J/kg*K
c_metal = (c_Cu * m_Cu + c_Al * m_Al)/(m_Cu + m_Al);
c_solid = 0.921e3; %J/kg*K  %silica gel

% Adsorption Model parameters
K = 2*10^-12; %Pa
Qst = 2.51 *10^6; %J/kg
% Qst = 2.693*10^6;
R = 8314; %J/kmol*K
MolarMass = meanMolarMass(water);
Rwater = R/MolarMass; %J/kg*K

%Pressures
setState_Tsat(water,[chilled_T_in,0]);
P_evap = pressure(water);
setState_Tsat(water,[cooling_T_in,1]);
P_cond = pressure(water);

% T1 = T_low;
% P1 = P_evap;
% set(water_vap,'P',P_evap,'T',T_low,'X','H2O:1');
%set(water_vap,'T',T1,'P',P1);
% density_amb = density(water_vap);

j = 1;
T_max = hot_T_in;

for cycle = 1:1
%-----------------------------------------------------------
    %HEATING
    T_evap = chilled_T_in;
    AdsorbA_Heating = Heating(hot_T_in);
    
%     if cycle > 1
%         TB_init = AdsorbB_Desorb.T_bed(end);
%         mB_init = AdsorbB_Desorb.m_gas(end);
%         mB_cond = AdsorbB_Desorb.m_cond(end);
%         AdsorbB_Cooling = Cooling(TB_init,mB_init,mB_cond);
%     end
%-----------------------------------------------------------
    %DESORPTION
    T_cond = cooling_T_in;  %condensation temperature
    TA_init = AdsorbA_Heating.T_bed(end);
    TA_final = hot_T_in;
    mA_init = AdsorbA_Heating.m_gas(end);
%     Thotstream = 79+273;%AdsorbA_Heating.heatingTout;
    AdsorbA_Desorb = Desorption(T_cond,TA_init,TA_final,mA_init);
    %     if cycle > 1
    %         TB_init = AdsorbB_Cooling.T_bed(end);
    %         TB_final = T_low;
    %         mB_init = AdsorbB_Cooling.m_gas(end);
    %         AdsorbB_Adsorb = Adsorption(TB_init,TB_final,mB_init,mB_cond);
    %     end
%-----------------------------------------------------------
    %CONDENSATION
    T2 = AdsorbA_Desorb.T_bed(end);
    m_ref = max(AdsorbA_Desorb.m_cond);
    Cond = Condenser(m_ref,T_max);
    Qcond = Cond.Qcond;
%-----------------------------------------------------------
    %Water states
    T1_water = AdsorbA_Heating.T_water(1) ;
    P1_water = P_evap;
    T2_water = AdsorbA_Heating.T_water(end);
    P2_water = P_cond;
    T3_water = T_cond;
    P3_water = P_cond;
    
    set(water,'T',T1_water,'P',P1_water);
    h1 = enthalpy_mass(water);
    set(water,'T',T2_water,'P',P2_water);
    h2 = enthalpy_mass(water);
    set(water,'T',T3_water,'P',P3_water);
    h3 = enthalpy_mass(water);
    h4 = h3;
    
%-----------------------------------------------------------
    %EXPANSION
    xin = Cond.xout;
    Exp = Expansion(h4,m_ref);
    T4_water = Exp.T;
    P4_water = P_evap;
%-----------------------------------------------------------
    %EVAPORATION
    % xin = Exp.xout;
    % y = Exp.y;
    Evap = Evaporator(Exp.hout,m_ref,T_low);
    Qevap = Evap.Qevap;
%-----------------------------------------------------------
    %COOLING
    TA_init = AdsorbA_Desorb.T_bed(end);
    mA_init = AdsorbA_Desorb.m_gas(end);
    mA_cond = AdsorbA_Desorb.m_cond(end);
    AdsorbA_Cooling = Cooling(TA_init,mA_init,mA_cond);
    
    %     AdsorbB_Heating = Heating(T_evap);
%-----------------------------------------------------------
    %ADSORPTION
    TA_init = AdsorbA_Cooling.T_bed(end);
    TA_final = T_low;
    mA_init = AdsorbA_Cooling.m_gas(end);
    Flowxin = Evap.Flowxout;
    AdsorbA_Adsorb = Adsorption(TA_init,TA_final,mA_init,mA_cond,Flowxin);
    
%     TB_init = AdsorbB_Heating.T_bed(end);
%     TB_final = hot_T_in;
%     mB_init = AdsorbB_Heating.m_gas(end);
%     AdsorbB_Desorb = Desorption(T_cond,TB_init,TB_final,mB_init);
%     
end

%-----------------------------------------------------------
%EXTERNAL PIPING STREAMS
[HotStream ChilStream CoolStream] = ExternalPiping(hot_T_in,chilled_T_in,cooling_T_in,hot_T_diff,chilled_T_diff,cooling_T_diff);
Qhotpaper = HotStream.Qheating;
mhotpaper = HotStream.m_dot;
hin = HotStream.h_in;
hout = HotStream.h_out;
Qhot = (AdsorbA_Heating.Q + AdsorbA_Desorb.Q_des)/cycletime;
mdot_hot = (-Qhot)/(hout-hin);
% fprintf('Qpaper = %.2f kW\n',Qhotpaper)
fprintf('Qhot = %.2f kW\n',Qhot/1e3)
% fprintf('mdotpaper = %.2f kg/s\n',mhotpaper)
% fprintf('mdot_hot = %.2f kg/s\n',mdot_hot)
% fprintf('***************************************************************\n')

Qchilpaper = ChilStream.Qchilled;
mchilpaper = ChilStream.m_dot;
hin = ChilStream.h_in;
hout = ChilStream.h_out;
Qchil = Evap.Qevap/cycletime;
mdot_chil = (-Qchil)/(hout-hin);
% fprintf('Qpaper = %.2f kW\n',Qchilpaper)
fprintf('Qchil = %.2f kW\n',Qchil/1e3)
% fprintf('mdotpaper = %.2f kg/s\n',mchilpaper)
% fprintf('mdot_chil = %.2f kg/s\n',mdot_chil)
% fprintf('***************************************************************\n')

Qcoolpaper = CoolStream.Qcooling;
mcoolpaper = CoolStream.m_dot;
hin = CoolStream.h_in;
hout = CoolStream.h_out;
Qcool = (AdsorbA_Cooling.Q + AdsorbA_Adsorb.Q + Qcond)/cycletime;
mdot_cool = (-Qcool)/(hout-hin);
% fprintf('Qpaper = %.2f kW\n',Qcoolpaper)
fprintf('Qcool = %.2f kW\n',Qcool/1e3)
% fprintf('mdotpaper = %.2f kg/s\n',mcoolpaper)
% fprintf('mdot_cool = %.2f kg/s\n',mdot_cool)
% fprintf('***************************************************************\n')

%-----------------------------------------------------------
%EXERGY CALCULATIONS
%Heating
HotFlowXin = mdot_hot*HotStream.FlowXin;
HotFlowXout = mdot_hot*HotStream.FlowXout;
DesFlowXout = AdsorbA_Desorb.FlowXout/cycletime;
DesDX = AdsorbA_Desorb.DX/cycletime;
XheatingDest = HotFlowXin - HotFlowXout - DesFlowXout - DesDX;

ChilFlowXin = mdot_chil*ChilStream.FlowXin;
ChilFlowXout = mdot_chil*ChilStream.FlowXout;
EvapFlowXin = Evap.FlowXin/(3*60);
EvapFlowXout = Evap.FlowXout/(3*60);
XchilDest = ChilFlowXin + EvapFlowXin - ChilFlowXout - EvapFlowXout;

CoolFlowXin = mdot_cool*CoolStream.FlowXin;
CoolFlowXout = mdot_cool*CoolStream.FlowXout;
AdsFlowXin = AdsorbA_Adsorb.FlowXin/cycletime;
AdsDX = AdsorbA_Adsorb.DX/cycletime;
XcoolDest = CoolFlowXin + AdsFlowXin - CoolFlowXout - AdsDX;

DXHot = HotFlowXin - HotFlowXout;
DXChil = ChilFlowXin - ChilFlowXout;
Xefficiency = -DXChil/DXHot;

% fprintf('Xefficiency = %.2f',Xefficiency)


[h_dome, s_dome, u_dome, T_dome, P_dome, v_dome] = vaporDome(water);

fprintf('\n***************************************************************\n')

Xheating = AdsorbA_Heating.Xloss/1e6;
Xdesorb = AdsorbA_Desorb.Xloss/1e6;
Xcond = Cond.Xloss/1e6;
Xcooling = AdsorbA_Cooling.Xloss/1e6;
Xevap = Evap.Xloss/1e6;
Xadsorb = AdsorbA_Adsorb.Xloss/1e6;
Xexpan = Exp.Xloss/1e6;
Xloss_system = Xheating + Xdesorb + Xcond + Xcooling + Xevap + Xadsorb + Xexpan;
Xefficiency = (Xevap)/(Xloss_system);

% fprintf('Xheating = %.2f MJ\n',Xheating)
% fprintf('Xdesorption = %.2f MJ\n',Xdesorb)
% fprintf('Xcondensation = %.2f MJ\n',Xcond)
% fprintf('Xcooling = %.2f MJ\n',Xcooling)
% fprintf('Xexpansion = %.2f MJ\n',Xexpan)
% fprintf('Xevaporation = %.2f MJ\n',Xevap)
% fprintf('Xadsorption = %.2f MJ\n',Xadsorb)
% fprintf('Xefficiency = %.2f',Xefficiency)

% fprintf('\n***************************************************************\n')

% Q_sens_evap = AdsorbA_Adsorb.Q_sens/1e6;
Q_heating = AdsorbA_Heating.Q/1e6;
Q_desorb = AdsorbA_Desorb.Q_des/1e6;
Q_cooling = AdsorbA_Cooling.Q/1e6;
Q_adsorb = AdsorbA_Adsorb.Q/1e6;
QevapMJ = Qevap/1e6;
QcondMJ = Qcond/1e6;
Qin = (Q_heating+Q_desorb);

% fprintf('Qheating = %.2f MJ\n',Q_heating)
% fprintf('Qdesorption = %.2f MJ\n',Q_desorb)
% fprintf('Qcondensation = %.2f MJ\n',QcondMJ)
% fprintf('Qcooling = %.2f MJ\n',Q_cooling)
% fprintf('Qevaporation = %.2f MJ\n',QevapMJ)
% fprintf('Qadsorption = %.2f MJ\n',Q_adsorb)

COP = Qchil/Qhot;
CoolingCap = (Qchil/1e3);
CoolingCap_perHot = Qevap/Qhot;
fprintf('COP = %.2f \n',COP)
% fprintf('Cooling Capacity = %.2f kW\n',CoolingCap)
% fprintf('Cooling Capacity = %.2f kJ/mass solid\n',CoolingCap_perHot)

Qadd = (Q_heating + Q_desorb + QevapMJ);
Qrej = -(Q_cooling + Q_adsorb + QcondMJ);
error = (Qadd - Qrej)/(Qadd) * 100;
% fprintf('error = %.2f percent\n',error)

% m_eff = m_solid * (q_max - q_min)
% m_util = 1 - (q_min/q_max)

%Figures
% Tticks = [10 20 30 40 50 60 70 80 90]+273';
% ticklabels = num2str(Tticks-273);
% Tlabels = {ticklabels};
% figure(1) % Ln P vs -1/T
% clf
% semilogy(-1./AdsorbA_Heating.T_bed,AdsorbA_Heating.P_bed/1e3,'bx-')
% hold on
% semilogy(-1./[T1_water T2_water T3_water T4_water T1_water],...
%     [P1_water P2_water P3_water P4_water P1_water]/1e3,'ro-')
% semilogy(-1./AdsorbA_Desorb.T_bed,AdsorbA_Desorb.P_bed/1e3,'bx-')
% semilogy(-1./AdsorbA_Cooling.T_bed,AdsorbA_Cooling.P_bed/1e3,'bx-')
% semilogy(-1./AdsorbA_Adsorb.T_bed,AdsorbA_Adsorb.P_bed/1e3,'bx-')
% hold off
% xlabel('Temperature (C)')
% ylabel('Pressure (kPa)')
% % legend('bed','water')
% set(gca,'XTick',-1./Tticks)
% set(gca,'XTickLabel',{Tticks-273})
% title('Duhring Plot for Adsorption Chiller Cycle')
% plotfixer
% 
% figure(2)
% clf
% hold on
% plot(AdsorbA_Heating.T_bed-273,AdsorbA_Heating.m_ads,'bx-')
% % plot(AdsorbB_Cooling.T_bed-273,AdsorbB_Cooling.m_gas,'rx-')
% plot(AdsorbA_Desorb.T_bed-273,AdsorbA_Desorb.m_ads,'gx-')
% % plot(AdsorbB_Adsorb.T_bed-273,AdsorbB_Adsorb.m_gas,'rx-')
% plot(AdsorbA_Cooling.T_bed-273,AdsorbA_Cooling.m_ads,'rx-')
% plot(AdsorbA_Adsorb.T_bed-273,AdsorbA_Adsorb.m_ads,'cx-')
% xlabel('T (C)')
% ylabel('Adsorbed water (kg)')
% legend('A','B')
% legend('Heat/Pressure','Desorb/Cond','Cool/Depressure','Evap/Adsorb')
% title('Adsorption Chiller Mass')
% text(25,225,'1')
% text(45,225,'2')
% text(85,25,'3')
% text(60,25,'4')
% plotfixer
% 
% 
% figure(3)
% clf
% hold on
% plot([h1 h2 h3 h4 h1]./1e3,...
%     [P1_water P2_water P3_water P4_water P1_water]./1e3,'rx-')
% plot(AdsorbA_Heating.h./1e3, AdsorbA_Heating.P_bed./1e3,'bx-')
% plot(AdsorbA_Desorb.h./1e3,AdsorbA_Desorb.P_bed./1e3,'bx-')
% plot(AdsorbA_Cooling.h./1e3,AdsorbA_Cooling.P_bed./1e3,'bx-')
% plot(AdsorbA_Adsorb.h./1e3, AdsorbA_Adsorb.P_bed./1e3,'bx-')
% text(h1/1e3,P1_water/1e3,'1')
% text(h2/1e3,P2_water/1e3,'2')
% text(h3/1e3,P3_water/1e3,'3')
% text(h4/1e3,P4_water/1e3,'4')
% plot(h_dome/1e3,P_dome./1e3,'k-')
% ylabel('P (kPa)')
% xlabel('Enthalpy (kJ/kg_{H2O} K)')
% title('Adsorption Chiller')
% legend('Refrigerant H2O','Vapor H2O in Bed')
% axis([-1.6e4 -1.3e4 1 5])
% plotfixer
% 
% figure(4) % Temperature vs Concentration
% clf
% hold on
% plot(AdsorbA_Heating.T_bed-273,AdsorbA_Heating.q,'x-')
% plot(AdsorbA_Desorb.T_bed-273,AdsorbA_Desorb.q,'x-')
% plot(AdsorbA_Cooling.T_bed-273,AdsorbA_Cooling.q,'x-')
% plot(AdsorbA_Adsorb.T_bed-273,AdsorbA_Adsorb.q,'x-')
% xlabel('T (C)')
% ylabel('Adsorbed water (kg_{H2O} / kg_{silica})')
% % legend('water','bed')
% title('Adsorption Chiller Concentrations')
% plotfixer

%Outputs
Outputs.COP = COP;
Outputs.CoolingCap = CoolingCap;
Outputs.error = error;

Outputs.Xheating = Xheating;
Outputs.Xdesorb = Xdesorb;
Outputs.Xcond = Xcond;
Outputs.Xcooling = Xcooling;
Outputs.Xevap = Xevap;
Outputs.Xadsorb = Xadsorb;
Outputs.Xexpan = Xexpan;

Outputs.Xsystem = Xloss_system;
Outputs.Xefficiency = Xefficiency;



end


