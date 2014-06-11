function[Outputs] = CrossOver(TempMode,EffMode,hot_T_in,chilled_T_in,cooling_T_in,hot_T_diff,chilled_T_diff,cooling_T_diff)

water = importPhase('liquidVapor.xml','water');

global open_volume porosity 
global m_solid m_metal c_metal c_solid
global P_evap P_cond T_max T_evap T_cond
global j cycletime
global To Po
global Qst K Rwater R

To = 300;
Po = 101325;

% cycletime = 6*60; %second %6minutes
t_normal=430; t_regen=30; t_htrec=20;
cycletime =(t_normal+t_regen+t_htrec);

evapEff = 1;
condEff = 1;
adsorberEff = 1;

%-----------------------------------------------------------
%Silica Gel Properties
%-----------------------------------------------------------
volume = 2.65; %3.1*.9*.9;  %m^3 chamber volume find value
adsorbent_density = 800; %800kg/m^3
porus_volume = 0.35; %ml/kg
water_content = 0; %percent
porosity = adsorbent_density*porus_volume*1000*1e-6; %fraction

%-----------------------------------------------------------
%Masses and Geometry
%-----------------------------------------------------------
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

%-----------------------------------------------------------
% Adsorption Model parameters
%-----------------------------------------------------------
K = 2*10^-12; %1/Pa
Qst = 2.51 *10^6; %J/kg
R = 8314; %J/kmol*K
MolarMass = meanMolarMass(water);
Rwater = R/MolarMass; %J/kg*K

%-----------------------------------------------------------
%Pressures
%-----------------------------------------------------------
setState_Tsat(water,[chilled_T_in,0]);
P_evap = pressure(water);
setState_Tsat(water,[cooling_T_in,1]);
P_cond = pressure(water);


j = 1;
T_max = hot_T_in;
T_evap = chilled_T_in;
T_cond = cooling_T_in;  %condensation temperature

switch EffMode
    case 0
        T_low = cooling_T_in;
    case 1
        T_low = cooling_T_in + cooling_T_diff;
end

%-----------------------------------------------------------
%CYCLE
%-----------------------------------------------------------
for cycle = 1:2
%-----------------------------------------------------------
%MODE A - Valves closed
%-----------------------------------------------------------
    %Adsorber A - Heating
    AdsorbA_Heating = Heating(T_low);
    AdA_Thotout = hot_T_in - hot_T_diff;
    
    %Adsorber B - Cooling
    if cycle > 1
        TB_init = AdsorbB_Desorb.T_bed(end);
        mgas_init = AdsorbB_Desorb.m_gas(end);
        AdsorbB_Cooling = Cooling(TB_init,mgas_init);
    end
    
%-----------------------------------------------------------
%MODE B - Valves opened
%-----------------------------------------------------------
    %Adsorber A - Desorption and Condensation
    Tdesorb = AdsorbA_Heating.T_bed(end);
    TA_final = hot_T_in;
    mgas_init = AdsorbA_Heating.m_gas(end);
    AdsorbA_Desorb = Desorption(T_cond,Tdesorb,TA_final,mgas_init,T_low);
    Amass_ref = AdsorbA_Desorb.m_cond;
    
    ACond = Condenser(Amass_ref,TA_final);
    AQcond = ACond.Qcond;
    
    %Adsorber B - Adsorption and Evaporation
    if cycle > 1
        BExp = Expansion(Bmass_ref);
        BEvap = Evaporator(Amass_ref,T_cond);
        BQevap = BEvap.Qevap; 
        
        TB_init = AdsorbB_Cooling.T_bed(end);
        TB_final = T_low;
        mgas_init = AdsorbB_Cooling.m_gas(end);
        Flowxin = BEvap.Flowxout;
        AdsorbB_Adsorb = Adsorption(TB_init,TB_final,mgas_init,Bmass_ref,Flowxin);
    end
   
%-----------------------------------------------------------
%MODE C - Valves closed
%-----------------------------------------------------------
    %Adsorber A - Cooling
    TA_init = AdsorbA_Desorb.T_bed(end);
    mgas_init = AdsorbA_Desorb.m_gas(end);
    AdsorbA_Cooling = Cooling(TA_init,mgas_init);
    
    %Adsorber B - Heating
    AdsorbB_Heating = Heating(T_low);
    AdB_Thotout = hot_T_in - hot_T_diff;
    
%-----------------------------------------------------------
%MODE D - Valves opened
%-----------------------------------------------------------
    %Adsorber A - Adsorption and Evaporation
    AExp = Expansion(Amass_ref);
    T4_water = AExp.T;
    AEvap = Evaporator(Amass_ref,T_cond);
    AQevap = AEvap.Qevap;
    
    Tadsorb = AdsorbA_Cooling.T_bed(end);
    TA_final = T_low;
    mgas_init = AdsorbA_Cooling.m_gas(end);
    Flowxin = AEvap.Flowxout;
    AdsorbA_Adsorb = Adsorption(Tadsorb,TA_final,mgas_init,Amass_ref,Flowxin);

    %Adsorber B - Desorption and Condensation
    TB_init = AdsorbB_Heating.T_bed(end);
    TB_final = hot_T_in;
    mgas_init = AdsorbB_Heating.m_gas(end);
    AdsorbB_Desorb = Desorption(T_cond,TB_init,TB_final,mgas_init,T_low);
    Bmass_ref = AdsorbA_Desorb.m_cond;
    
    BCond = Condenser(Bmass_ref,TA_final);
    BQcond = BCond.Qcond;
    
%-----------------------------------------------------------  
%Water states
%-----------------------------------------------------------
    T1_water = AdsorbA_Heating.T_water(1) ;
    P1_water = P_evap;
    T2_water = AdsorbA_Heating.T_water(end);
    P2_water = P_cond;
    T3_water = T_cond;
    P3_water = P_cond;
    P4_water = P_evap;
    
    set(water,'T',T1_water,'P',P1_water);
    h1 = enthalpy_mass(water);
    set(water,'T',T2_water,'P',P2_water);
    h2 = enthalpy_mass(water);
    set(water,'T',T3_water,'P',P3_water);
    h3 = enthalpy_mass(water);
    h4 = h3;
   
end

%-----------------------------------------------------------
%EXTERNAL PIPING STREAMS
%-----------------------------------------------------------
[HotStream ChilStream CoolStream] = ExternalPiping(hot_T_in,chilled_T_in,cooling_T_in,hot_T_diff,chilled_T_diff,cooling_T_diff);
Qhotpaper = HotStream.Qheating;
mhotpaper = HotStream.m_dot;
hin = HotStream.h_in;
hout = HotStream.h_out;
Qh = AdsorbA_Heating.Q;
Qdes = AdsorbA_Desorb.Q_des;
Qdothot = adsorberEff*(AdsorbA_Heating.Q + AdsorbA_Desorb.Q_des)/cycletime;
mdot_hot = (-Qdothot)/(hout-hin);
% fprintf('Qpaper = %.2f kW\n',Qhotpaper)
% fprintf('Qhot = %.2f kW\n',Qdothot/1e3)
% fprintf('mdotpaper = %.2f kg/s\n',mhotpaper)
% fprintf('mdot_hot = %.2f kg/s\n',mdot_hot)
% fprintf('***************************************************************\n')

Qchilpaper = ChilStream.Qchilled;
mchilpaper = ChilStream.m_dot;
hin = ChilStream.h_in;
hout = ChilStream.h_out;
Dh = hout-hin;
Qdotchil = AQevap/cycletime;
mdot_chil = (-Qdotchil)/(hout-hin);
% fprintf('Qpaper = %.2f kW\n',Qchilpaper)
% fprintf('Qchil = %.2f kW\n',Qdotchil/1e3)
% fprintf('mdotpaper = %.2f kg/s\n',mchilpaper)
% fprintf('mdot_chil = %.2f kg/s\n',mdot_chil)
% fprintf('***************************************************************\n')

Qcoolpaper = CoolStream.Qcooling;
mcoolpaper = CoolStream.m_dot;
hin = CoolStream.h_in;
hout = CoolStream.h_out;
Qdotcool = (AdsorbA_Cooling.Q + AdsorbA_Adsorb.Q)/cycletime;
mdot_cool = (-Qdotcool)/(hout-hin);

mdot_cond = (AQcond/cycletime)/(hout-hin);
% fprintf('Qpaper = %.2f kW\n',Qcoolpaper)
% fprintf('Qcool = %.2f kW\n',Qdotcool/1e3)
% fprintf('mdotpaper = %.2f kg/s\n',mcoolpaper)
% fprintf('mdot_cool = %.2f kg/s\n',mdot_cool)
% fprintf('***************************************************************\n')

%-----------------------------------------------------------
%EXERGY CALCULATIONS
%-----------------------------------------------------------
Xviolation = 0;

HotFlowXin = mdot_hot*HotStream.Flowxin;
HotFlowXout = mdot_hot*HotStream.Flowxout;
DesFlowXout = AdsorbA_Desorb.FlowXout/cycletime;
DesDX = AdsorbA_Desorb.DX/cycletime;
XheatingDest = HotFlowXin - HotFlowXout - DesFlowXout - DesDX;
if XheatingDest < 0
    fprintf('!!Heating exergy destruction violation. Tempurature out of range \n')
    Xviolation = 1;
end

ChilFlowXin = mdot_chil*ChilStream.Flowxin;
ChilFlowXout = mdot_chil*ChilStream.Flowxout;
EvapFlowXin = AEvap.FlowXin/(6*60);
EvapFlowXout = AEvap.FlowXout/(6*60);
XchilDest = ChilFlowXin + EvapFlowXin - ChilFlowXout - EvapFlowXout;
if XchilDest < 0
    fprintf('!!Chilling exergy destruction violation. Tempurature out of range \n')
    Xviolation = 1;
end

CoolFlowXin = mdot_cool*CoolStream.Flowxin;
CoolFlowXout = mdot_cool*CoolStream.Flowxout;
DXCoolFlow = CoolFlowXin - CoolFlowXout;
AdsFlowXin = AdsorbA_Adsorb.FlowXin/cycletime;
AdsDX = AdsorbA_Adsorb.DX/cycletime;

CondFlowXin = ACond.FlowXin/cycletime;
CondFlowXout = ACond.FlowXout/cycletime;

% XcoolDest = CoolFlowXin + CondFlowXin + AdsFlowXin...
%     - CoolFlowXout - CondFlowXout - AdsDX

XcoolDest = CoolFlowXin + AdsFlowXin - CoolFlowXout - AdsDX;
if XcoolDest < 0
    fprintf('!!Cooling from adsorber exergy destruction violation. Tempurature out of range \n')
    Xviolation = 1;
end

CoolCondFlowXin = mdot_cond*CoolStream.Flowxin;
CoolCondFlowXout = mdot_cond*CoolStream.Flowxout;
XcondDest = CoolCondFlowXin + CondFlowXin - CoolCondFlowXout - CondFlowXout;
if XcondDest < 0
    fprintf('!!Cooling from condenser exergy destruction violation. Tempurature out of range \n')
    Xviolation = 1;
end

XexpanDest = AExp.XexpDest/cycletime;
if XexpanDest < 0
    fprintf('!!Flash exergy destruction violation.\n')
    Xviolation = 1;
    AExp.y
end

DXHot = HotFlowXin - HotFlowXout;
DXChil = ChilFlowXin - ChilFlowXout;
Xefficiency_HotReuse = -DXChil/DXHot;
Xefficiency = -DXChil/HotFlowXin;
fprintf('Exergy efficiency with hot reuse = %.2f\n',Xefficiency_HotReuse)
fprintf('Exergy efficiency = %.2f\n',Xefficiency)

XsystemSum = XheatingDest + XchilDest + XcoolDest + XexpanDest;
Xsystem = HotFlowXin + ChilFlowXin + CoolFlowXin... 
    - HotFlowXout - ChilFlowXout - CoolFlowXout;

DiffX = Xsystem - XsystemSum;
Xerror = DiffX/Xsystem;
fprintf('System exergy error = %.2f\n',Xerror)

%Other Chiller Comparisons
COPcompress = 5;
Wcompress = Qdotchil/COPcompress;
Xcomp = Wcompress + EvapFlowXout - CondFlowXin;
XeffCompression = -DXChil/Xcomp;

% fprintf('\n***************************************************************\n')

%-----------------------------------------------------------
%ENERGY CALCULATIONS
%-----------------------------------------------------------
QheatingMJ = AdsorbA_Heating.Q/1e6;
QdesorbMJ = AdsorbA_Desorb.Q_des/1e6;
QcoolingMJ = AdsorbA_Cooling.Q/1e6;
QadsorbMJ = AdsorbA_Adsorb.Q/1e6;
QevapMJ = AQevap/1e6;
QcondMJ = AQcond/1e6;

% fprintf('Qheating = %.2f MJ\n',QheatingMJ)
% fprintf('Qdesorption = %.2f MJ\n',QdesorbMj)
% fprintf('Qcondensation = %.2f MJ\n',QcondMJ)
% fprintf('Qcooling = %.2f MJ\n',QcoolingMJ)
% fprintf('Qevaporation = %.2f MJ\n',QevapMJ)
% fprintf('Qadsorption = %.2f MJ\n',QadsorbMJ)

% Qadd = (Q_heating + Q_desorb + QevapMJ);
% Qrej = -(Q_cooling + Q_adsorb + QcondMJ);
% error = (Qadd - Qrej)/(Qadd) * 100;
% fprintf('error = %.2f percent\n',error)

%-----------------------------------------------------------
%PERFORMANCE METRICS
%-----------------------------------------------------------
Tads_avg = (Tadsorb + T_low)/2;
Tdes_avg = (T_max + Tdesorb)/2;
Tevap_avg = (T_low + T_evap)/2;

COP = Qdotchil/Qdothot;
CoolingCap = (Qdotchil/1e3);
HeatingLoad = Qdothot/1e3;
CoolingCap_perHot = CoolingCap/mdot_hot;
Cooling_perHot = mdot_chil/mdot_hot;
HCR = (m_solid*c_solid)/(m_metal*c_metal);
SCP = (AQevap/1e3)/m_solid;
COP_carnot = (-(Tads_avg - Tdes_avg)/Tdes_avg)*(Tevap_avg/(Tads_avg-Tevap_avg));

fprintf('COP = %.2f \n',COP)
fprintf('Cooling Capacity = %.2f kW\n',CoolingCap)
fprintf('Heating Load = %.2f kW\n',HeatingLoad)
% fprintf('Cooling Capacity = %.2f mdot_chil/mdot_hot\n',Cooling_perHot)
% fprintf('Cooling Capacity = %.2f kW Chil/(kg/s)Hot\n',CoolingCap_perHot)
fprintf('Heat Capacity Ratio (silica/metal)= %.2f \n',HCR)
fprintf('Specific Cooling Power = %.2f kJ/kg_silica\n',SCP)
fprintf('COP carnot = %.2f \n',COP_carnot)
fprintf('Efficiency = %.2f \n',COP/COP_carnot)


%-----------------------------------------------------------
% FIGURES
%-----------------------------------------------------------
if TempMode == 1
    Tticks = [10 20 30 40 50 60 70 80 90]+273';
    ticklabels = num2str(Tticks-273);
    Pticks = [1 2 3 4 5];
    figure() % Ln P vs -1/T
    clf
    semilogy(-1./AdsorbA_Heating.T_bed,AdsorbA_Heating.P_bed/1e3,'bx-')
    hold on
    semilogy(-1./[T1_water T2_water T3_water T4_water T1_water],...
        [P1_water P2_water P3_water P4_water P1_water]/1e3,'ro-')
    semilogy(-1./AdsorbA_Desorb.T_bed,AdsorbA_Desorb.P_bed/1e3,'bx-')
    semilogy(-1./AdsorbA_Cooling.T_bed,AdsorbA_Cooling.P_bed/1e3,'bx-')
    semilogy(-1./AdsorbA_Adsorb.T_bed,AdsorbA_Adsorb.P_bed/1e3,'bx-')
    hold off
    xlabel('Temperature (\circC)')
    ylabel('Pressure (kPa)')
    legend('Adsorption bed','Refrigerant water','Location','SouthEast')
    set(gca,'XTick',-1./Tticks)
    set(gca,'XTickLabel',{Tticks-273})
    set(gca,'YTick',Pticks)
    set(gca,'YTickLabel',{Pticks})
    axis([-1/283 -1/373 1 5])
    text(-1/(30+273),1.4,'1')
    text(-1/(50+273),4.7,'2')
    text(-1/(87+273),4.7,'3')
    text(-1/(65+273),1.6,'4')
    text(-1/(31+273),4.7,'3w')
    text(-1/(12+273),1.4,'4w')
    title('Duhring Plot for Adsorption Chiller Cycle')
    plotfixer
    
    figure()
    clf
    hold on
    plot(AdsorbA_Heating.T_bed-273,AdsorbA_Heating.m_ads,'bx-')
    % plot(AdsorbB_Cooling.T_bed-273,AdsorbB_Cooling.m_gas,'rx-')
    plot(AdsorbA_Desorb.T_bed-273,AdsorbA_Desorb.m_ads,'gx-')
    % plot(AdsorbB_Adsorb.T_bed-273,AdsorbB_Adsorb.m_gas,'rx-')
    plot(AdsorbA_Cooling.T_bed-273,AdsorbA_Cooling.m_ads,'rx-')
    plot(AdsorbA_Adsorb.T_bed-273,AdsorbA_Adsorb.m_ads,'cx-')
    xlabel('T (\circC)')
    ylabel('Adsorbed water (kg)')
    legend('A','B')
    legend('Heat/Pressure','Desorb/Cond','Cool/Depressure','Evap/Adsorb')
    title('Adsorption Chiller Mass')
    text(31,170,'1')
    text(50,170,'2')
    text(85,25,'3')
    text(60,25,'4')
    plotfixer
    
    [h_dome, s_dome, u_dome, T_dome, P_dome, v_dome] = vaporDome(water);
    figure()
    clf
    hold on
    plot([h1 h2 h3 h4 h1]./1e3 - h1/1e3,...
        [P1_water P2_water P3_water P4_water P1_water]./1e3,'rx-')
    plot(AdsorbA_Heating.h./1e3 - h1/1e3, AdsorbA_Heating.P_bed./1e3,'bx-')
    plot(AdsorbA_Desorb.h./1e3 - h1/1e3,AdsorbA_Desorb.P_bed./1e3,'bx-')
    plot(AdsorbA_Cooling.h./1e3- h1/1e3,AdsorbA_Cooling.P_bed./1e3,'bx-')
    plot(AdsorbA_Adsorb.h./1e3- h1/1e3, AdsorbA_Adsorb.P_bed./1e3,'bx-')
    text(h1/1e3- h1/1e3,P1_water/1e3,'1')
    text(h2/1e3- h1/1e3,P2_water/1e3,'2')
    text(h3/1e3- h1/1e3,P3_water/1e3,'3w')
    text(h4/1e3- h1/1e3,P4_water/1e3,'4w')
    plot(h_dome/1e3- h1/1e3,P_dome./1e3,'k-')
    ylabel('Pressure (kPa)')
    xlabel('Enthalpy (kJ/kg_{H2O}) relative to state 1')
    % title('Adsorption Chiller')
    legend('Refrigerant H_2O','Vapor H_2O in Bed','location','Best')
    axis([-2500 200 1 5])
    plotfixer
    
    figure() % Temperature vs Concentration
    clf
    hold on
    plot(AdsorbA_Heating.T_bed-273,AdsorbA_Heating.q,'bx-')
    plot(AdsorbA_Desorb.T_bed-273,AdsorbA_Desorb.q,'gx-')
    plot(AdsorbA_Cooling.T_bed-273,AdsorbA_Cooling.q,'rx-')
    plot(AdsorbA_Adsorb.T_bed-273,AdsorbA_Adsorb.q,'cx-')
    xlabel('T (\circC)')
    ylabel('Adsorbed water (kg_{H2O} / kg_{silica})')
    legend('Heat/Pressure','Desorb/Cond','Cool/Depressure','Evap/Adsorb')
    title('Adsorption Chiller Concentrations')
    text(30,.19,'1')
    text(50,.19,'2')
    text(85,.03,'3')
    text(60,.03,'4')
    plotfixer
end

%-----------------------------------------------------------
%OUTPUTS
%-----------------------------------------------------------
Outputs.COP = COP;
Outputs.CoolingCap = CoolingCap;
% Outputs.error = error;

Outputs.Xchill = XchilDest;
Outputs.Xcooling = XcoolDest;
Outputs.Xheating = XheatingDest;
Outputs.Xexpan = XexpanDest;
Outputs.Xsystem = XsystemSum;
Outputs.Xefficiency = Xefficiency_HotReuse;
Outputs.Xviolation = Xviolation;

Outputs.Qdothot = Qdothot;
Outputs.Qdotchil = Qdotchil;

Outputs.qmax = AdsorbA_Heating.q;
Outputs.qmin = AdsorbA_Cooling.q;
Outputs.Pheating = AdsorbA_Heating.P_bed;

Outputs.mr = Amass_ref;


end


