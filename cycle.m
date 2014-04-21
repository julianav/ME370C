clear all;

water = importPhase('liquidVapor.xml','water');
water_vap = IdealGasMix('gri30.xml');

T_amb = 300;

%Properties
volume = 1*.75*1.5;  %m^3 chamber volume find value
adsorbent_density = 800; %kg/m^3
porus_volume = 0.35; %ml/kg
water_content = 0; %percent
porosity = adsorbent_density*porus_volume*1000*1e-6; %fraction

m_metal = 0;
c_metal = 0;   %metal cover
c_solid = 0.921; %kJ/kg*K  %silica gel


% Adsorption Model parameters
K = 5.5*10^-12; %Pa
Qst = 2.37 *10^3;
R = 8.314; %J/mol*K
MolarMass = meanMolarMass(water);

%Hot water stream
hot_V_dot = 18; %m^3/h
density_water = 968;
hot_m_dot = hot_V_dot*density_water/3600;
hot_T_in = 85 + 273; %K
hot_T_diff = 5.6; %K
hot_T_out = hot_T_in + hot_T_diff;
hot_P_in = 10e5; %Pa
hot_P_drop = 44100; %Pa = 4.5 mH20
hot_P_out = hot_P_in - hot_P_drop;
set(water,'P',hot_P_in,'T',hot_T_in);
hot_h_in = enthalpy_mass(water);
set(water,'P',hot_P_out,'T',hot_T_out);
hot_h_out = enthalpy_mass(water);

%Chilled water stream 
chilled_V_dot = 12; %m^3/h
chilled_m_dot = chilled_V_dot*density_water/3600;
chilled_T_in = 14 + 273; %K
chilled_T_out = 9 + 273; %K
chilled_P_in = 10e5; %Pa 
chilled_P_drop = 53900; %Pa = 5.5 mH20
chilled_P_out = chilled_P_in - chilled_P_drop;
set(water,'P',chilled_P_in,'T',chilled_T_in);
chilled_h_in = enthalpy_mass(water);
set(water,'P',chilled_P_out,'T',chilled_T_in);
chilled_h_out = enthalpy_mass(water);

%Cooling water stream
cooling_V_dot = 42;
cooling_m_dot = cooling_V_dot*density_water/3600;
cooling_T_in = 31 + 273; %K
cooling_T_diff = 3.8;
cooling_T_out = cooling_T_in + cooling_T_diff;
cooling_P_in = 10e5; %Pa
cooling_P_drop = 78453.2; %Pa = 8 mH20
cooling_P_out = cooling_P_in - cooling_P_drop;
set(water,'P',cooling_P_in,'T',cooling_T_in);
cooling_h_in = enthalpy_mass(water);
set(water,'P',cooling_P_out,'T',cooling_T_out);
cooling_h_out = enthalpy_mass(water);

%Pressures
setState_Tsat(water,[chilled_T_out,0]);
P_evap = pressure(water);
setState_Tsat(water,[cooling_T_in,1]);
P_cond = pressure(water);

% P_evap = 1.2e3; %Pa %approximations
% P_cond = 4.5e3; %Pa

%Masses
m_solid = volume * adsorbent_density; %kg    % silica_Mass = 895; %kg

T1 = T_amb;
P1 = P_evap;
set(water_vap,'P',P_evap,'T',T_amb,'X','H2O:1');
%set(water_vap,'T',T1,'P',P1);
density_amb = density(water_vap);
% m_total = volume * density_amb * porosity;

% skeletal_density = 0.1785; %kg/m^3  helium density
% Ep = 1 - (adsorbent_density/skeletal_density);  %partical porosity
% porosity = Ep + (1 - Ep) * porus_volume;
% total_porus_volume = chamber_volume*adsorb_density*porus_volume*1000*1e-6; %m^3
% 
% m_adsorbed = total_porus_volume*water_content;
% m_vapor = total_porus_volume/density_water;
% 
% m_total = m_adsorbed + m_vapor;



%-----------------------------------------------------
%Process 1 = Heating and Pressurization
%-----------------------------------------------------
T_evap = chilled_T_in;
setState_Tsat(water,[T_evap 1]);
P_evap = pressure(water);

% set(water,'T',T_amb,'P',P_evap);
set(water_vap,'P',P_evap,'T',T_amb,'X','H2O:1');
ug1 = intEnergy_mass(water_vap);
hg1 = enthalpy_mass(water_vap);
rhog1 = density(water_vap);

setState_Tsat(water,[T_amb 0]);
rhoL = density(water);
rho_a1 = rhoL;

q_max = K*exp(Qst/(R/MolarMass*T_amb))*P_evap;

P1_water = P_evap;
T1_water = T_amb;
P1_bed = P_evap;
T1_bed = T_amb;

% P2_water = P_cond;
% T2_water = Qst/(R/MolarMass*log(q_max/(K*P_cond)));
% P2_bed = P_cond;
% T2_bed = T2_water;

dP = (P_cond - P_evap)/10;
i = 1;
T_prev = T1_water;
ug_prev = ug1;
ha = hg1 - Qst;  %enthalpy of adsorbate phase
rho_a = rhoL;
ua_prev = ha - P_evap / rho_a;
Q12 = 0;

for P = P_evap:dP:P_cond
    T = Qst/(R/MolarMass*log(q_max/(K*P)));
    
    set(water_vap,'P',P,'T',T,'X','H2O:1');
    hg = enthalpy_mass(water_vap);
    rhog = density(water_vap);
    ug = intEnergy_mass(water_vap);
    
    setState_Tsat(water,[T 0]);
    rhoL = density(water);
    
    %---------------------------------------------
    dT = T - T_prev;
    dQ_solid = m_solid * c_solid * dT;
    %---------------------------------------------
    dtheta = 0;
    dW = volume * P * dtheta;
    
    m_a = m_solid * q_max;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P / rho_a;
    du_a = u_a - ua_prev;
    
    dQ_adsorbate = dW * m_a * du_a;
    %----------------------------------------------
%     m_g = m_total - m_a;
    m_g = (volume * density_amb) * (porosity - (q_max * m_solid)/(volume * rho_a1));
    du_g = ug - ug_prev;
    dQ_gas = -dW + m_g * du_g;
    %-----------------------------------------------
    dQ12 = dQ_solid + dQ_adsorbate + dQ_gas;
    %-----------------------------------------------
    Q12 = Q12 + dQ12;
    
    %graphing
    T12_water(i) = T;
    P12_water(i) = P;
    q12(i) = q_max;
    
    ua_prev = u_a;
    T_prev = T;
    i = i + 1;
end

T12_bed = T12_water;
P12_bed = P12_water;

% % setState_Psat(water,[P_cond 0]);
% % T2 = temperature(water);
% set(water,'T',T2_water,'P',P_cond);
% hg2 = enthalpy_mass(water);
% % vg2 = 1/density(water);
% ug2 = intEnergy_mass(water);
% 
% dT = T2_water - T1_water;
% dQ_solid = m_solid * c_solid * dT;
% dtheta = 0;
% dW = volume * P1 * dtheta;
% m_a = m_solid * q_max;
% du_a = ug2 - ug1;
% du_g = du_a;
% dQ_adsorbate = dW * m_a * du_a;
% 
% m_g = m_total - m_a;
% dQ_gas = -dW + m_g * du_g;
% 
% Q12 = dQ_solid + dQ_adsorbate + dQ_gas

Q_heating = hot_m_dot*(hot_h_out - hot_h_in);  %reference value only


%%--------------------------------------------------------------
%Process 2 = Desorption and Condensation;  %constant P
%%--------------------------------------------------------------
T_cond = cooling_T_in;  %condensation temperature
setState_Tsat(water,[T_cond 1]); %Vapor phase
P_cond = pressure(water);  %condensation pressure
hg2 = enthalpy_mass(water);
ug2 = intEnergy_mass(water);

T_cond = cooling_T_in;  %condensation temperature
setState_Tsat(water,[T_cond 0]); %liquid phase
rhoL = 995; %density(water2);

T2_bed = T12_bed(end);
T3_water = T_cond;
P3_water = P_cond;
T3_bed = hot_T_in;
P3_bed = P_cond;

%q_min = K*exp(Qst/(R/MolarMass*T_cond))*P_cond;

dT = (T3_bed - T2_bed)/10;
q_prev = q_max;
ha = hg2 - Qst;  %enthalpy of adsorbate phase
rho_a = rhoL;
ua_prev = ha - P_cond / rho_a;
ug_prev = ug2;
Q_cond = 0;
Q23 = 0;
m_des = 0;
i = 1;

for t = T2_bed:dT:T3_bed
    T = t;
    set(water,'T',T,'P',P_cond);
    hg = enthalpy_mass(water);
    ug = intEnergy_mass(water); 
    %-----------------------------------------
    dQ_solid = m_solid * c_solid * dT;
    %-----------------------------------------
    q_curr = K*exp(Qst/(R/MolarMass*T))*P_cond;
    dq = q_curr - q_prev;
    m_a = q_curr * m_solid;
    dm_a = dq * m_solid;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P_cond / rho_a;
    du_a = u_a - ua_prev;
    dU_a = m_a * du_a + u_a * dm_a;
    
    dm_desorbed = -dm_a;
    
    dW_a = volume * P_cond * dq;
    
    dQ_adsorbate = dU_a + hg * dm_desorbed + dW_a;
    %---------------------------------------------
    %     m_g = m_total - m_a;
    m_g = (volume * density_amb) * (porosity - (q_max * m_solid)/(volume * rho_a));
    dm_g = -dm_desorbed;
    du_g = ug - ug_prev;
    
    dU_g = m_g * du_g + ug * dm_g;
    
    dQ_gas = dU_g - dW_a - dm_g * hg;
    %---------------------------------------------
    dQ23 = dQ_solid + dQ_adsorbate + dQ_gas;
    %---------------------------------------------
    dm_cond = -(dm_a + dm_g);
    dQ_cond = hg * dm_cond;
    %---------------------------------------------
    Q23 = Q23 + dQ23;
    Q_cond = Q_cond + dQ_cond;
    
    m_des = m_des + dm_desorbed;
    
    %graphing
    T23_water(i) = T;  %?????
    P23_water(i) = P_cond;
    T23_bed(i) = T;
    P23_bed(i) = P_cond;
    q23(i) = q_curr;
    m_desorbed(i) = m_des;

    ua_prev = u_a;
    ug_prev = ug;
    q_prev = q_curr;
    
    i = i + 1;
end
 
% h_out = cooling_h_in + Q_cond / cooling_m_dot;
% set(water,'P',P_cond,'H',h_out);
T_cooling_out = temperature(water);

Q_cooling = cooling_m_dot*(cooling_h_out - cooling_h_in);



%State 2
% m_dot_water = cooling_m_dot;
% h2 = h3-Q_cond/m_dot_water;
% P2 = P_cond;
% set(water,'P',P2,'H',h2);
% T2 = temperature(water);

%Process 3 = Cooling and Depressurization
T_max = T3_bed;
P4_bed = P_evap;

set(water_vap,'P',P_evap,'T',T_amb,'X','H2O:1');
ug3 = intEnergy_mass(water_vap);
hg3 = enthalpy_mass(water_vap);
rhog3 = density(water_vap);

setState_Psat(water,[P_evap 0]);
T4_bed = temperature(water);
rhoL = density(water);
rho_a3 = rhoL;

q_min = K*exp(Qst/(R/MolarMass*T_max))*P_cond;

dP = -(P_cond - P_evap)/10;
T_prev = T_max;
ug_prev = ug3;
ha = hg3 - Qst;  %enthalpy of adsorbate phase
rho_a = rhoL;
ua_prev = ha - P_evap / rho_a;
i = 1;
Q34 = 0;

for P = P_cond:dP:P_evap
    T = Qst/(R/MolarMass*log(q_min/(K*P)));
    
    set(water_vap,'P',P,'T',T,'X','H2O:1');
    hg = enthalpy_mass(water_vap);
    rhog = density(water_vap);
    ug = intEnergy_mass(water_vap);
    
    setState_Tsat(water,[T 0]);
    rhoL = density(water);
    
    %---------------------------------------------------
    dT = T - T_prev;
    dQ_metal = -m_metal * c_metal * dT;
    dQ_solid = -m_solid * c_solid * dT;
    %---------------------------------------------------
    dtheta = 0; %Va/Vb;
    m_a = m_solid * q_min;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P / rho_a;
    du_a = u_a - ua_prev;
    
    dQ_adsorbate = -volume * P * dtheta - m_a * du_a;
    %---------------------------------------------------
    m_g = (volume * rhog3) * (porosity - (q_min * m_solid)/(volume * rho_a3));
    du_g = ug - ug_prev;
    
    dQ_gas = volume * P * dtheta - m_g * du_g;
    %---------------------------------------------------
    dQ34 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
    %---------------------------------------------------
    Q34 = Q34 + dQ34;
    
    %graphing
    T34_bed(i) = T;
    P34_bed(i) = P;
    q34(i) = q_min;
    
    T_prev = T;
    ua_prev = u_a;
    ug_prev = ug;
    i = i + 1;
end




%Process 4




% [h_dome, s_dome, u_dome, T_dome, P_dome, v_dome] = VaporDome(water);


%% Figures

figure(1)
clf
hold on
% plot([T1_water T2_water T3_water],[P1_water P2_water P3_water],'x-')
% plot([T1_bed T2_bed T3_bed],[P1_bed P2_bed P3_bed ],'x-')

plot(T12_bed-273,P12_bed/1e3,'x-')
plot(T23_bed-273,P23_bed/1e3,'x-')
plot(T34_bed-273,P34_bed/1e3,'x-')
xlabel('T(C)')
ylabel('P(kPa)')
title('Adsorption Chiller Cycle')
% legend('water','bed')
% axis([200 400 1000 5000])
hold off

% figure(2)
% clf
% semilogy([-1/T1_water -1/T2_water -1/T3_water],[P1_water P2_water P3_water],'bx-')
% hold on
% semilogy([-1/T1_bed -1/T2_bed -1/T3_bed],[P1_bed P2_bed P3_bed],'ro-')
% axis([200 400 1000 5000])
% hold off
% xlabel('-1/T')
% ylabel('lnP')
% legend('water','bed')
% title('Adsorption Chiller Cycle')

figure(3)
clf
hold on
plot(T12_bed-273,q12,'x-')
plot(T23_bed-273,q23,'x-')
plot(T34_bed-273,q34,'x-')
xlabel('T (C)')
ylabel('q (kg/kg)')
% legend('water','bed')
title('Adsorption Chiller Concentrations')

figure(4)
clf
hold on
plot(T23_bed-273,m_desorbed,'x-')
xlabel('T (C)')
ylabel('mass desorbed (kg)')
% legend('water','bed')
title('Adsorption Chiller Mass Desorbed')
axis([45 90 0 180])
