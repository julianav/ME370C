clear all
format compact
fprintf('\n***************************************************************\n')
fprintf('NEW RUN')
fprintf('\n***************************************************************\n')

water = importPhase('liquidVapor.xml','water');
water_vap = IdealGasMix('gri30.xml');

global open_volume porosity 
global m_solid m_metal c_metal c_solid Qst
global P_evap P_cond T_amb T_max
global j
global To Po
global dead
global Qst K

T_amb = 300;

To = 300;
Po = 101325;
% set(water_vap,'T',dead.To,'P',dead.Po,'X','H2O');
% set(water,'T',dead.To,'P',dead.Po);
% dead.mu_o = chemPotentials(water);
% dead.so = entropy_mass(water);
% dead.uo = intEnergy_mass(water);
% dead.vo = 1/density(water);

%Properties
volume = 3.3*.9*.9;  %m^3 chamber volume find value
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
R = 8314; %J/kmol*K
MolarMass = meanMolarMass(water);
Rwater = R/MolarMass; %J/kg*K

%Hot water stream
hot_T_in = 85 + 273; %K

%Chilled water stream 
chilled_T_in = 14 + 273; %K
chilled_T_out = 9 + 273; %K

%Cooling water stream
cooling_T_in = 31 + 273; %K
cooling_V_dot = 42;
density_water = 968;
mdot_cooling = cooling_V_dot*density_water/3600;

%Pressures
setState_Tsat(water,[chilled_T_out,0]);
P_evap = pressure(water);
setState_Tsat(water,[cooling_T_in,1]);
P_cond = pressure(water);

T1 = T_amb;
P1 = P_evap;
set(water_vap,'P',P_evap,'T',T_amb,'X','H2O:1');
%set(water_vap,'T',T1,'P',P1);
density_amb = density(water_vap);

j = 1;
T_max = hot_T_in;

for cycle = 1:1
    cycle
    %Step 1
    T_evap = chilled_T_in;
    AdsorbA_Heating = Heating(T_evap);
    if cycle > 1
%         TB_init = AdsorbB_Desorb.T_bed(end);
%         mB_init = AdsorbB_Desorb.m_gas(end);
%         mB_cond = AdsorbB_Desorb.m_cond(end);
%         AdsorbB_Cooling = Cooling(TB_init,mB_init,mB_cond);
    end
    
    %Step 2
    T_cond = cooling_T_in;  %condensation temperature
    TA_init = AdsorbA_Heating.T_bed(end);
    TA_final = hot_T_in;
    mA_init = AdsorbA_Heating.m_gas(end);
    AdsorbA_Desorb = Desorption(T_cond,TA_init,TA_final,mA_init);
    
    if cycle > 1
%         TB_init = AdsorbB_Cooling.T_bed(end);
%         TB_final = T_amb;
%         mB_init = AdsorbB_Cooling.m_gas(end);
%         AdsorbB_Adsorb = Adsorption(TB_init,TB_final,mB_init,mB_cond);
    end
    
    %Step 3
    TA_init = AdsorbA_Desorb.T_bed(end);
    mA_init = AdsorbA_Desorb.m_gas(end);
    mA_cond = AdsorbA_Desorb.m_cond(end);
    AdsorbA_Cooling = Cooling(TA_init,mA_init,mA_cond);
    
%     AdsorbB_Heating = Heating(T_evap);
    
    %Step 4
    TA_init = AdsorbA_Cooling.T_bed(end);
    TA_final = T_amb;
    mA_init = AdsorbA_Cooling.m_gas(end);
    AdsorbA_Adsorb = Adsorption(TA_init,TA_final,mA_init,mA_cond);
%     
%     TB_init = AdsorbB_Heating.T_bed(end);
%     TB_final = hot_T_in;
%     mB_init = AdsorbB_Heating.m_gas(end);
%     AdsorbB_Desorb = Desorption(T_cond,TB_init,TB_final,mB_init);
    
end

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

m_ref = AdsorbA_Adsorb.m_evap(end);

%Condensation
xin = AdsorbA_Desorb.xin;
Cond = Condenser(m_ref,mdot_cooling,cooling_T_in,xin);
Qcond2 = Cond.Qcond/1e6


%Throttling
xin = Cond.xout;
Exp = Expansion( h4,m_ref,xin );
T4_water = Exp.T;
P4_water = P_evap;

%Evaporation
xin = Exp.xout;
xout = AdsorbA_Adsorb.xin; 
Evap = Evaporator(Exp.hfg,m_ref,xin,xout,cooling_T_in);
Qevap = Evap.Qevap;


[h_dome, s_dome, u_dome, T_dome, P_dome, v_dome] = vaporDome(water);

fprintf('\n***************************************************************\n')

Xdesorb = AdsorbA_Desorb.Xloss
Xcond = Cond.Xloss
Xevap = Evap.Xloss
Xadsorb = AdsorbA_Adsorb.Xloss

Xloss_system = Xdesorb + Xcond + Xevap + Xadsorb

Q_sens_evap = AdsorbA_Adsorb.Q_sens/1e6

Q_heating = AdsorbA_Heating.Q/1e6
Q_desorb = AdsorbA_Desorb.Q_des/1e6
Q_cooling = AdsorbA_Cooling.Q/1e6
Q_adsorb = AdsorbA_Adsorb.Q/1e6
Qevap = Qevap/1e6;
Qcond = AdsorbA_Desorb.Q_cond/1e6
Qin = (Q_heating+Q_desorb)

COP = Qevap/Qin

Qadd = (Q_heating + Q_desorb + Qevap);
Qrej = -(Q_cooling + Q_adsorb - Qcond2);
error = abs((Qadd - Qrej)/Qadd * 100);
fprintf('error = %.2f percent\n',error)

% m_eff = m_solid * (q_max - q_min)
% m_util = 1 - (q_min/q_max)

%%
%Figures
Tticks = [10 20 30 40 50 60 70 80 90]+273';
ticklabels = num2str(Tticks-273);
Tlabels = {ticklabels};
figure(1) % Ln P vs -1/T
clf
semilogy(-1./AdsorbA_Heating.T_bed,AdsorbA_Heating.P_bed,'bx-')
hold on
semilogy(-1./[T1_water T2_water T3_water T4_water T1_water],...
    [P1_water P2_water P3_water P4_water P1_water],'ro-')
semilogy(-1./AdsorbA_Desorb.T_bed,AdsorbA_Desorb.P_bed,'bx-')
semilogy(-1./AdsorbA_Cooling.T_bed,AdsorbA_Cooling.P_bed,'bx-')
semilogy(-1./AdsorbA_Adsorb.T_bed,AdsorbA_Adsorb.P_bed,'bx-')
hold off
xlabel('Temperature (K)')
ylabel('Pressure (MPa)')
% legend('bed','water')
set(gca,'XTick',-1./Tticks)
set(gca,'XTickLabel',{Tticks-273})
title('Duhring Plot for Adsorption Chiller Cycle')
plotfixer

figure(2)
clf
hold on
plot(AdsorbA_Heating.T_bed-273,AdsorbA_Heating.m_ads,'bx-')
% plot(AdsorbB_Cooling.T_bed-273,AdsorbB_Cooling.m_gas,'rx-') 
plot(AdsorbA_Desorb.T_bed-273,AdsorbA_Desorb.m_ads,'gx-')
% plot(AdsorbB_Adsorb.T_bed-273,AdsorbB_Adsorb.m_gas,'rx-')
plot(AdsorbA_Cooling.T_bed-273,AdsorbA_Cooling.m_ads,'rx-')
plot(AdsorbA_Adsorb.T_bed-273,AdsorbA_Adsorb.m_ads,'cx-')
xlabel('T (C)')
ylabel('mass (kg)')
legend('A','B')
legend('Heat/Pressure','Desorb/Cond','Cool/Depressure','Evap/Adsorb')
title('Adsorption Chiller Mass')
text(25,225,'1')
text(45,225,'2')
text(85,25,'3')
text(60,25,'4')
plotfixer


figure(3)
clf
hold on
plot([h1 h2 h3 h4 h1]./1e3,...
    [P1_water P2_water P3_water P4_water P1_water]./1e3,'rx-')
plot(AdsorbA_Heating.h./1e3, AdsorbA_Heating.P_bed./1e3,'bx-')
plot(AdsorbA_Desorb.h./1e3,AdsorbA_Desorb.P_bed./1e3,'bx-')
plot(AdsorbA_Cooling.h./1e3,AdsorbA_Cooling.P_bed./1e3,'bx-')
plot(AdsorbA_Adsorb.h./1e3, AdsorbA_Adsorb.P_bed./1e3,'bx-')
text(h1/1e3,P1_water/1e3,'1')
text(h2/1e3,P2_water/1e3,'2')
text(h3/1e3,P3_water/1e3,'3')
text(h4/1e3,P4_water/1e3,'4')
plot(h_dome/1e3,P_dome./1e3,'k-')
ylabel('P (kPa)')
xlabel('Enthalpy (kJ/kg K)')
title('Adsorption Chiller')
legend('Refrigerant','Bed')
axis([-1.6e4 -1.3e4 1 5])
plotfixer

figure(4) % Temperature vs Concentration
clf
hold on
plot(AdsorbA_Heating.T_bed-273,AdsorbA_Heating.q,'x-')
plot(AdsorbA_Desorb.T_bed-273,AdsorbA_Desorb.q,'x-')
plot(AdsorbA_Cooling.T_bed-273,AdsorbA_Cooling.q,'x-')
plot(AdsorbA_Adsorb.T_bed-273,AdsorbA_Adsorb.q,'x-')
xlabel('T (C)')
ylabel('q (kg/kg)')
% legend('water','bed')
title('Adsorption Chiller Concentrations')
plotfixer





%%
% %-----------------------------------------------------
% %Process 1-2 = Heating and Pressurization 
% %-----------------------------------------------------
% 
% 
% setState_Tsat(water,[T_evap 1]);  %vapor
% P_evap = pressure(water);
% 
% set(water,'T',T_amb,'P',P_evap);
% ug1 = intEnergy_mass(water);
% hg1 = enthalpy_mass(water);
% rhog1 = density(water);
% mu1 = chemPotentials(water);
% X1 = exergy_mass(water,To,Po,mu_o);
% 
% setState_Tsat(water,[T_amb 0]);
% rhoL = density(water);
% rho_a1 = rhoL;
% 
% q_max = Adsorbate_Con_Ratio(T_amb,P_evap);
% 
% m_a = q_max * m_solid;
% m_gas1 = open_volume * rhog1 * porosity;
% m_total = m_a + m_gas1;
% 
% P1_water = P_evap;
% T1_water = T_amb;
% P1_bed = P_evap;
% T1_bed = T_amb;
% 
% dP = (P_cond - P_evap)/10;
% i = 1;
% j = 1;
% T_prev = T1_water;
% ug_prev = ug1;
% ha = hg1 - Qst;  %enthalpy of adsorbate phase
% rho_a = rhoL;
% ua_prev = ha - P_evap / rho_a;
% Q12 = 0;
% 
% for P = P_evap:dP:P_cond
%     T = Qst/(Rwater*log(q_max/(K*P)));
%     
%     set(water,'P',P,'T',T);
% %     set(water_vap,'P',P,'T',T,'X','H2O:1');
%     hg = enthalpy_mass(water);
%     rhog = density(water);
%     ug = intEnergy_mass(water);
%     
%     setState_Tsat(water,[T 0]);  %liquid
%     rhoL = density(water);
%     
%     %---------------------------------------------
%     dT = T - T_prev;
%     dQ_metal = m_metal * c_metal * dT;
%     dQ_solid = m_solid * c_solid * dT;
%     %---------------------------------------------
%     dtheta = 0;
%     dW = - volume * P * dtheta;
%     
%     m_a = q_max * m_solid;
%     
%     ha = hg - Qst;  %enthalpy of adsorbate phase
%     rho_a = rhoL;
%     u_a = ha - P / rho_a;
%     du_a = u_a - ua_prev;
%     
%     dQ_adsorbate = m_a * du_a - dW;
%     %----------------------------------------------
%     m_g = m_total - m_a;
% %     m_g_paper = (volume * density_amb) * (porosity - (q_max * m_solid)/(volume * rho_a1));
%     du_g = ug - ug_prev;
%     dQ_gas = -dW + m_g * du_g;
%     %-----------------------------------------------
%     dQ12 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
%     %-----------------------------------------------
%     Q12 = Q12 + dQ12;
%     
%     %graphing
%     T12_water(i) = T;
%     P12_water(i) = P;
%     T12_bed(i) = T;
%     P12_bed(i) = P;
%     q12(i) = q_max;
%     m_gas12(i) = m_g;
%     m_ads12(i) = m_a;
%     Q12_heat(i) = Q12;
%     h12(i) = hg;
%     m_ref(j) = 0;
%     P_ref(j) = P;
%     
%     ua_prev = u_a;
%     T_prev = T;
%     i = i + 1;
%     j = j + 1;
% end
% 
% 
% % % setState_Psat(water,[P_cond 0]);
% % % T2 = temperature(water);
% % set(water,'T',T2_water,'P',P_cond);
% % hg2 = enthalpy_mass(water);
% % % vg2 = 1/density(water);
% % ug2 = intEnergy_mass(water);
% % 
% % dT = T2_water - T1_water;
% % dQ_solid = m_solid * c_solid * dT;
% % dtheta = 0;
% % dW = volume * P1 * dtheta;
% % m_a = m_solid * q_max;
% % du_a = ug2 - ug1;
% % du_g = du_a;
% % dQ_adsorbate = dW * m_a * du_a;
% % 
% % m_g = m_total - m_a;
% % dQ_gas = -dW + m_g * du_g;
% % 
% % Q12 = dQ_solid + dQ_adsorbate + dQ_gas
% 
% % Q_heating_ref = hot_m_dot*(hot_h_out - hot_h_in)/1e6  %reference value only
% 
% 
% %%--------------------------------------------------------------
% %Process 2-3 = Desorption and Condensation;  %constant P
% %%--------------------------------------------------------------
% % T_cond = cooling_T_in;  %condensation temperature
% % setState_Tsat(water,[T_cond 1]); %Vapor phase
% % set(water,'T',T_cond,'P',P_cond);
% % % P_cond = pressure(water);  %condensation pressure
% % hg2 = enthalpy_mass(water);
% % ug2 = intEnergy_mass(water);
% 
% T_cond = cooling_T_in;  %condensation temperature
% setState_Tsat(water,[T_cond 0]); %liquid phase
% rhoL = density(water);
% 
% setState_Tsat(water,[T_cond,0]);
% hf_low = enthalpy_mass(water);
% 
% T2_bed = T12_bed(end);
% T2_water = T2_bed;
% P2_water = P_cond;
% P2_bed = P_cond;
% 
% T3_water = T_cond;
% P3_water = P_cond;
% T3_bed = hot_T_in;
% P3_bed = P_cond;
% 
% % q_min = K*exp(Qst/(R/MolarMass*T_cond))*P_cond;
% m_a2 = q_max * m_solid;
% m_g2 = m_gas12(end);
% 
% %Initial state
% set(water,'T',T2_bed,'P',P_cond);
% % P_cond = pressure(water);  %condensation pressure
% hg2 = enthalpy_mass(water);
% ug2 = intEnergy_mass(water);
% 
% 
% dT = (T3_bed - T2_bed)/10;
% q_prev = q_max;
% ha = hg2 - Qst;  %enthalpy of adsorbate phase
% rho_a = rhoL;
% ua_prev = ha - P_cond / rho_a;
% ug_prev = ug2;
% ma_prev = m_a2;
% mg_prev = m_g2;
% T_prev = T2_bed;
% Qcond = 0;
% Q23 = 0;
% Qgas = 0;
% m_des = 0;
% m_leaving_prev = 0;
% m_cond = 0;
% i = 1;
% 
% for T = T2_bed:dT:T3_bed
%     set(water,'T',T,'P',P_cond);
% %     set(water_vap,'P',P_cond,'T',T,'X','H2O:1');
%     hg = enthalpy_mass(water);
%     ug = intEnergy_mass(water);
%     rhog = density(water);
%     
%     mg_capacity = open_volume * rhog * porosity;
% 
%     setState_Tsat(water,[T 0]); %liquid
%     rhoL = density(water);
% 
%     %-----------------------------------------
%     dT = T - T_prev;
%     dQ_metal = m_metal * c_metal * dT;
%     dQ_solid = m_solid * c_solid * dT;
%     %-----------------------------------------
%     q = K*exp(Qst/(Rwater*T))*P_cond;
%     Dq = q - q_prev;
%     
%     m_a = q * m_solid;
%     Dm_a = m_a - ma_prev;
%     
%     ha = hg - Qst;  %enthalpy of adsorbate phase
%     rho_a = rhoL;
%     u_a = ha - P_cond / rho_a;
%     du_a = u_a - ua_prev;
%     dU_a = m_a * du_a + u_a * Dm_a;
%     
%     Dm_des = -Dm_a;
%     
%     dW_a = - volume * P_cond * Dq;
%     
%     dQ_adsorbate = dU_a - dW_a - Dm_a * hg;  %changed mass of des to mass of ads
%     %---------------------------------------------
% %     m_g_paper = (volume * rhog) * (porosity - (q * m_solid)/(volume * rho_a));
%     Dm_g = Dm_des;
% %     m_g_t = m_total - m_a;
% 
%     mg_init = mg_prev + Dm_g; % mass of gas initially
%     
%     if mg_init > mg_capacity;   % excess mass leaves the adsorber
%         m_leaving = mg_init - mg_capacity;
%     else
%         m_leaving = 0;
%     end
%     
%     m_g = mg_init - m_leaving;
%     Dm_g = m_g - mg_prev;
%     
%     du_g = ug - ug_prev;
%     dU_g = m_g * du_g + ug * Dm_g;  %m_g = Dm_a - m_leaving
%     
%     dW_g = -dW_a;
%     
%     Dm_cond = m_leaving;
%     m_cond = m_cond + Dm_cond;
% %     Dm_cond = Dm_g;
% %     dq_cond = Dm_cond * (hg);
%     dQ_cond = Dm_cond * (hg - hf_low);
%     
%     
%     dQ_gas = dU_g - dW_g - Dm_g * hg; %+ dQ_cond;
%     Qgas = Qgas + dQ_gas;
%     %---------------------------------------------
%     dQ23 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
%     %---------------------------------------------
% 
%     %---------------------------------------------
%     Q23 = Q23 + dQ23;
%     Qcond = Qcond + dQ_cond;
%     
%     m_des = m_des + Dm_des;
%     
%     %graphing
%     T23_water(i) = T;  %?????
%     P23_water(i) = P_cond;
%     T23_bed(i) = T;
%     P23_bed(i) = P_cond;
%     q23(i) = q;
%     m_desorbed(i) = m_des;
%     m_gas23(i) = m_g;
%     m_ads23(i) = m_a;
%     m_cond23(i) = m_cond;
%     Q23_des(i) = Q23;
%     Q23_gas(i) = Qgas;
%     h23(i) = hg;
%     m_ref(j) = m_cond;
%     P_ref(j) = P_cond;
% 
%     T_prev = T;
%     ua_prev = u_a;
%     ug_prev = ug;
%     q_prev = q;
%     ma_prev = m_a;
%     mg_prev = m_g;
%     m_leaving_prev = m_leaving;
%     
%     i = i + 1;
%     j = j + 1;
% end
%  
% % h_out = cooling_h_in + Q_cond / cooling_m_dot;
% % set(water,'P',P_cond,'H',h_out);
% T_cooling_out = temperature(water);
% 
% % Q_cooling_ref = cooling_m_dot*(cooling_h_out - cooling_h_in)/1e6
% 
% 
% 
% %State 2
% % m_dot_water = cooling_m_dot;
% % h2 = h3-Q_cond/m_dot_water;
% % P2 = P_cond;
% % set(water,'P',P2,'H',h2);
% % T2 = temperature(water);
% 
% %---------------------------------------------------
% %Process 3-4 = Cooling and Depressurization
% %----------------------------------------------------
% T_max = T3_bed;
% P4_bed = P_evap;
% 
% % set(water_vap,'P',P_evap,'T',T_amb,'X','H2O:1');
% set(water,'T',T_max,'P',P_cond);
% ug3 = intEnergy_mass(water);
% hg3 = enthalpy_mass(water);
% rhog3 = density(water);
% 
% setState_Psat(water,[P_evap 0]);  %liquid
% T4_bed = temperature(water);
% rhoL = density(water);
% rho_a3 = rhoL;
% 
% q_min = K*exp(Qst/(Rwater*T_max))*P_cond;
% 
% dP = -(P_cond - P_evap)/10;
% T_prev = T_max;
% ug_prev = ug3;
% ha = hg3 - Qst;  %enthalpy of adsorbate phase
% rho_a = rhoL;
% ua_prev = ha - P_evap / rho_a;
% i = 1;
% Q34 = 0;
% 
% for P = P_cond:dP:P_evap
%     T = Qst/(Rwater*log(q_min/(K*P)));
%     
%     set(water_vap,'P',P,'T',T,'X','H2O:1');
%     hg = enthalpy_mass(water_vap);
%     rhog = density(water_vap);
%     ug = intEnergy_mass(water_vap);
%     
%     setState_Tsat(water,[T 0]);
%     rhoL = density(water);
%     
%     %---------------------------------------------------
%     dT = T - T_prev;
%     dQ_metal = m_metal * c_metal * dT;
%     dQ_solid = m_solid * c_solid * dT;
%     %---------------------------------------------------
%     dtheta = 0; %Va/Vb;
%     dW = - volume * P * dtheta;
%     
%     m_a = q_min * m_solid;
%     
%     ha = hg - Qst;  %enthalpy of adsorbate phase
%     rho_a = rhoL;
%     u_a = ha - P / rho_a;
%     du_a = u_a - ua_prev;
%     
% %     dQ_adsorbate = -volume * P * dtheta - m_a * du_a;
%     dQ_adsorbate = m_a * du_a - dW;
%     %---------------------------------------------------
%     m_g_paper = (volume * rhog3) * (porosity - (q_min * m_solid)/(volume * rho_a3));
%     m_g = m_gas23(end);
%     du_g = ug - ug_prev;
%     
%     dQ_gas = m_g * du_g - dW;
%     %---------------------------------------------------
%     dQ34 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
%     %---------------------------------------------------
%     Q34 = Q34 + dQ34;
%     
%     %graphing
%     T34_bed(i) = T;
%     P34_bed(i) = P;
%     q34(i) = q_min;
%     m_gas34(i) = m_g;
%     m_ads34(i) = m_a;
%     Q34_cool(i) = Q34;
%     h34(i) = hg;
%     m_ref(j) = m_cond;
%     P_ref(j) = P;
%     
%     T_prev = T;
%     ua_prev = u_a;
%     ug_prev = ug;
%     i = i + 1;
%     j = j + 1;
% end
% 
% 
% %-------------------------------------------------------------
% %Process 4-1
% %-------------------------------------------------------------
% 
% T4_bed = T34_bed(end);
% set(water_vap,'P',P_evap,'T',T4_bed,'X','H2O:1');
% ug4 = intEnergy_mass(water_vap);
% hg4 = enthalpy_mass(water_vap);
% rhog4 = density(water_vap);
% 
% setState_Tsat(water,[T4_bed 0]); %liquid phase
% rhoL = density(water);
% 
% m_a4 = q_min * m_solid;
% m_g4 = m_gas34(end);
% 
% ha = hg4 - Qst;  %enthalpy of adsorbate phase
% rho_a = rhoL;
% ua_prev = ha - P_cond / rho_a;
% q_prev = q_min;
% ug_prev = ug4;
% ma_prev = m_a4;
% mg_prev = m_g4;
% T_prev = T4_bed;
% Q41 = 0;
% m_ads = 0;
% m_evap = 0;
% Q_evap = 0;
% i = 1;
% dT = -(T4_bed - T1_bed)/10;
% 
% for T = T4_bed:dT:T1_bed
%     q = K*exp(Qst/(Rwater*T))*P_evap;
%     
%     set(water,'T',T,'P',P_evap);
% %     set(water_vap,'P',P_evap,'T',T,'X','H2O:1');
%     hg = enthalpy_mass(water);
%     rhog = density(water);
%     ug = intEnergy_mass(water);
% 
%     mg_capacity = open_volume * rhog * porosity;
%     
%     setState_Psat(water,[P_evap,1]);
%     h_ev = enthalpy_mass(water);
%     
%     setState_Tsat(water,[T 0]);
%     rhoL = density(water);
%     
%     %--------------------------------------------------
% %     dT = T - T_prev;
%     dQ_metal = m_metal * c_metal * dT;
%     dQ_solid = m_solid * c_solid * dT;
%     %--------------------------------------------------
%     Dq = q - q_prev;
%     m_a = q * m_solid;
%     Dm_a = m_a - ma_prev;
%     
%     ha = hg - Qst;  %enthalpy of adsorbate phase
%     rho_a = rhoL;
%     u_a = ha - P_evap / rho_a;
%     du_a = u_a - ua_prev;
%     dU_a = m_a * du_a + u_a * Dm_a;
% 
%     Dm_ads = Dm_a;
%     
%     dW_a = - volume * P_evap * Dq;
%     
%     dQ_adsorbate = -h_ev * Dm_ads + dU_a - dW_a;
%     %-------------------------------------------------
% %     m_g_paper = (volume * density_amb) * (porosity - (q_max * m_solid)/(volume * rho_a1));
% 
%    
%     m_g = mg_capacity;
%     Dm_g = m_g - mg_prev;
%     Dm_entering = Dm_ads - Dm_g;
%  
%     du_g = ug - ug_prev;
%     
%     dU_g = m_g * du_g + ug * Dm_g;
%     
%     dQ_gas = dU_g - dW_a - Dm_g * h_ev;
%     %--------------------------------------------------
%     dQ41 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
%     %--------------------------------------------------
%     Dm_evap = Dm_entering;
%     m_evap = m_evap + Dm_evap;
%     
%     dQ_evap = Dm_evap * (h_ev);
%     %--------------------------------------------------
%     Q41 = Q41 + dQ41;
%     Q_evap = Q_evap + dQ_evap;
%     
%     m_ads = m_ads + Dm_ads;
%     
%     %graphing
% %     T41_water(i) = T;  %?????
%     P41_water(i) = P_evap;
%     T41_bed(i) = T;
%     P41_bed(i) = P_evap;
%     q41(i) = q;
%     m_adsorbed(i) = m_ads;
%     m_gas41(i) = m_g;
%     m_ads41(i) = m_a;
%     m_evap41(i) = m_evap;
%     Q41_ads(i) = Q41;
%     h41(i) = hg;
%     m_ref(j) = m_cond - m_evap;
%     P_ref(j) = P_evap;
%   
%     T_prev = T;
%     q_prev = q;
%     ua_prev = u_a;
%     ug_prev = ug;
%     ma_prev = m_a;
%     mg_prev = m_g;
%     
%     i = i + 1;
%     j = j + 1;
% end
% 
% %Condensation
% 
% % % Throttling
% % P4_water = P_evap;
% % setState_Psat(water,[P_evap,0]); %liquid
% % hf = enthalpy_mass(water);
% % T_evap = temperature(water);
% % T4_water = T_evap;
% % 
% % setState_Psat(water,[P_evap,1]); %vapor
% % hg = enthalpy_mass(water);
% % 
% % hfg = hf - hg;
% % 
% % set(water,'T',T_amb,'P',P_evap);
% % hf_amb = enthalpy_mass(water);
% % hf_evap = hf;
% % 
% % y = 1 - (hf_amb - hf_evap) / hfg;
% % 
% % %Evaporation
% % Qevap = m_evap * hfg;
% 
% %Water states
% set(water,'T',T1_water,'P',P1_water);
% h1 = enthalpy_mass(water);
% set(water,'T',T2_water,'P',P2_water);
% h2 = enthalpy_mass(water);
% set(water,'T',T3_water,'P',P3_water);
% h3 = enthalpy_mass(water);
% h4 = h3;
% 
% %Throttling
% setState_Psat(water,[P_evap,0]); %liquid
% hf = enthalpy_mass(water);
% setState_Psat(water,[P_evap,1]); %vapor
% hg = enthalpy_mass(water);
% T4_water = temperature(water);
% P4_water = P_evap;
% x = (h4 - hf)/(hg - hf);
% y = 1-x;
% 
% %Evaporation
% hfg = hg - hf;
% Qevap = y * m_evap * (hg - hf);
% 
% 
% [h_dome, s_dome, u_dome, T_dome, P_dome, v_dome] = vaporDome(water);
% 
% fprintf('\n***************************************************************\n')
% m_cond
% m_evap
% 
% Q12 = Q12/1e6
% Q23 = Q23/1e6
% Q34 = Q34/1e6
% Q41 = Q41/1e6
% Qevap = Qevap/1e6
% Qcond = Qcond/1e6
% Qin = (Q12+Q23)
% 
% COP = Qevap/Qin
% 
% Qadd = (Q12 + Q23 + Qevap)
% Qrej = -(Q34 + Q41 - Qcond)
% 
% 
% %% Figures
% 
% figure(1) % Temperature vs Pressure
% clf
% hold on
% plot([T12_water T3_water T4_water T1_water]-273,...
%     [P12_water P3_water P4_water P1_water]/1e3,'co-')
% plot(chilled_T_in-273,P_evap/1e3,'ro')
% plot(T12_bed-273,P12_bed/1e3,'x-')
% plot(T23_bed-273,P23_bed/1e3,'x-')
% plot(T34_bed-273,P34_bed/1e3,'x-')
% plot(T41_bed-273,P41_bed/1e3,'x-')
% xlabel('T(C)')
% ylabel('P(kPa)')
% title('Adsorption Chiller Cycle')
% % legend('water','bed')
% % axis([200 400 1000 5000])
% hold off
% 
% 
% Tticks = [10 20 30 40 50 60 70 80 90]';
% ticklabels = num2str(Tticks);
% Tlabels = {ticklabels};
% figure(2) % Ln P vs -1/T
% clf
% semilogy(-1./T12_bed,P12_bed,'bx-')
% hold on
% semilogy(-1./[T12_water T3_water T4_water T1_water],...
%     [P12_water P3_water P4_water P1_water],'ro-')
% semilogy(-1./T23_bed,P23_bed,'bx-')
% semilogy(-1./T34_bed,P34_bed,'bx-')
% semilogy(-1./T41_bed,P41_bed,'bx-')
% % axis([200 400 1000 5000])
% hold off
% xlabel('Temperature (K)')
% ylabel('Pressure (MPa)')
% legend('bed','water')
% set(gca,'XTick',-1 ./ Tticks)
% set(gca,'XTickLabel',Tlabels)
% title('Duhring Plot for Adsorption Chiller Cycle')
% 
% % figure(3) % Temperature vs Concentration
% % clf
% % hold on
% % plot(T12_bed-273,q12,'x-')
% % plot(T23_bed-273,q23,'x-')
% % plot(T34_bed-273,q34,'x-')
% % plot(T41_bed-273,q41,'x-')
% % xlabel('T (C)')
% % ylabel('q (kg/kg)')
% % % legend('water','bed')
% % title('Adsorption Chiller Concentrations')
% 
% % figure(4)
% % clf
% % hold on
% % plot(T12_bed-273,m_ads12,'bx-')
% % plot(T23_bed-273,m_ads23,'rx-')
% % plot(T34_bed-273,m_ads34,'gx-')
% % plot(T41_bed-273,m_ads41,'kx-')
% % % plot(T12_bed-273,m_gas12,'bo-')
% % % % plot([T12_bed T23_bed]-273,[m_ads12 + m_gas12 m_ads23 + m_gas23],'g-')
% % % plot(T23_bed-273,m_gas23,'ro-')
% % % plot(T34_bed-273,m_gas34,'go-')
% % % plot(T41_bed-273,m_gas41,'ko-')
% % xlabel('T (C)')
% % ylabel('mass (kg)')
% % % legend('adsorbed phase','gas phase','total')
% % legend('Heat/Pressure','Desorb/Cond','Cool/Depressure','Evap/Adsorb')
% % title('Adsorption Chiller Mass')
% 
% figure(5)
% clf
% hold on
% plot(T12_bed-273,m_gas12,'bo-')
% plot(T23_bed-273,m_gas23,'ro-')
% plot(T34_bed-273,m_gas34,'go-')
% plot([T34_bed(end) T41_bed(1)]-273,[m_gas34(end) m_gas41(1)],'ko-')
% plot(T41_bed-273,m_gas41,'ko-')
% xlabel('T (C)')
% ylabel('Vapor Mass (kg)')
% legend('Heat/Pressure','Desorb/Cond','Cool/Depressure','Evap/Adsorb')
% title('Vapor Mass in Adsorber')
% 
% 
% 
% % figure(6)
% % clf
% % hold on
% % plot(T12_bed-273,Q12_heat/1e6,'bx-')
% % plot(T23_bed-273,Q23_des/1e6,'bx-')
% % xlabel('T (C)')
% % ylabel('Driving energy (MJ)')
% % title('Adsorption Chiller')
% % % axis([45 90 0 180])
% % 
% % figure(7)
% % clf
% % hold on
% % plot(T34_bed-273,-Q34_cool/1e6,'bx-')
% % plot(T41_bed-273,-Q41_ads/1e6,'bx-')
% % xlabel('T (C)')
% % ylabel('Driving energy (MJ)')
% % title('Adsorption Chiller')
% % % axis([45 90 0 180])
% 
% % figure(8)
% % clf
% % hold on
% % plot(T23_bed-273,Q23_des/1e6,'bx-')
% % plot(T23_bed-273,Q23_gas/1e6,'rx-')
% % xlabel('T (C)')
% % ylabel('Driving energy (MJ)')
% % title('Adsorption Chiller')
% % legend('Total','Desorbed')
% % % axis([45 90 0 180])
% 
% figure(9)
% clf
% hold on
% plot([h1 h2 h3 h4 h1]./1e3,...
%     [P1_water P2_water P3_water P4_water P1_water]./1e3,'rx-')
% plot(h12./1e3,P12_bed./1e3,'bx-')
% plot(h23./1e3,P23_bed./1e3,'bx-')
% plot(h34./1e3,P34_bed./1e3,'bx-')
% plot(h41./1e3,P41_bed./1e3,'bx-')
% text(h1/1e3,P1_water/1e3,'1')
% text(h2/1e3,P2_water/1e3,'2')
% text(h3/1e3,P3_water/1e3,'3')
% text(h4/1e3,P4_water/1e3,'4')
% plot(h_dome/1e3,P_dome./1e3,'k-')
% ylabel('P (kPa)')
% xlabel('Enthalpy (kJ/kg K)')
% title('Adsorption Chiller')
% legend('Refrigerant','Bed')
% axis([-1.6e4 -1.3e4 1 5])
% 
% % 
% % figure(10)
% % clf
% % hold on
% % plot(T23_bed-273,m_cond23,'bo-')
% % plot(T41_bed-273,m_evap41,'go-')
% % xlabel('Temperature of Bed (C)')
% % ylabel('Refrigerant Mass (kg)')
% % legend('Mass Condensed','Mass Evaporated')
% % text(T2_bed-273,m_cond23(1),'2')
% % text(T3_bed-273,m_cond23(end),'3')
% % text(T4_bed-273,m_evap41(1),'4')
% % text(T1_bed-273,m_evap41(end),'1')
% % title('Refrigerant Mass')
% 
% figure(10)
% clf
% hold on
% plot(P_ref/1e3,m_ref,'bo-')
% % plot(P23_bed/1e3,m_cond23,'bo-')
% % plot(P41_bed/1e3,m_evap41,'go-')
% xlabel('Pressure Bed (kPa)')
% ylabel('Refrigerant Mass (kg)')
% text(P2_bed/1e3,m_cond23(1),'2')
% text(P3_bed/1e3,m_cond23(end),'3')
% text(P4_bed/1e3,m_evap41(1),'1')
% text(P1_bed/1e3,m_evap41(end),'4')
% text(3.75,80,'Condensation')
% text(3,140,'Throttling')
% text(1.7,80,'Evaporation')
% text(3,10,'Heating')
% title('Refrigerant Mass')
% 
% plotfixer