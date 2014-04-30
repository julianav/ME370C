clear all
format compact
fprintf('\n***************************************************************\n')
fprintf('NEW RUN')
fprintf('\n***************************************************************\n')

water = importPhase('liquidVapor.xml','water');
water_vap = IdealGasMix('gri30.xml');

T_amb = 300;

%Properties
volume = 3.4*1*1.5;  %m^3 chamber volume find value
adsorbent_density = 800; %800kg/m^3
porus_volume = 0.35; %ml/kg
water_content = 0; %percent
porosity = adsorbent_density*porus_volume*1000*1e-6; %fraction

%Masses
Ntubes = 672;
SA_tubes = 3.4 * (pi * 0.02) * Ntubes;
m_metal = SA_tubes * .8e-3 * 8940; %volume * density 
volume_tube = pi * .02^2 * 3.4 * Ntubes;
open_volume = volume - volume_tube;
m_solid = (open_volume) * adsorbent_density - m_metal; %kg    % silica_Mass = 895; %kg
c_metal = 0.39e3; %J/kg*K  %metal tubes, copper
c_solid = 0.921e3; %J/kg*K  %silica gel

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
setState_Tsat(water,[chilled_T_in,0]);
P_evap = pressure(water);
setState_Tsat(water,[cooling_T_in,1]);
P_cond = pressure(water);

% P_evap = 1.2e3; %Pa %approximations
% P_cond = 4.5e3; %Pa

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
%Process 1-2 = Heating and Pressurization 
%-----------------------------------------------------
T_evap = chilled_T_in;
setState_Tsat(water,[T_evap 1]);  %vapor
P_evap = pressure(water);

T_cool = cooling_T_in;
set(water,'T',T_amb,'P',P_evap);
% set(water_vap,'P',P_evap,'T',T_amb,'X','H2O:1');
ug1 = intEnergy_mass(water);
hg1 = enthalpy_mass(water);
rhog1 = density(water);

setState_Tsat(water,[T_amb 0]);
rhoL = density(water);
rho_a1 = rhoL;

q_max = K*exp(Qst/(R/MolarMass*T_amb))*P_evap;  %mass adsorbate / mass solid

m_a = q_max * m_solid;
m_gas1 = open_volume * rhog1 * porosity;
m_total = m_a + m_gas1;


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
j = 1;
T_prev = T1_water;
ug_prev = ug1;
ha = hg1 - Qst;  %enthalpy of adsorbate phase
rho_a = rhoL;
ua_prev = ha - P_evap / rho_a;
Q12 = 0;

for P = P_evap:dP:P_cond
    T = Qst/(R/MolarMass*log(q_max/(K*P)));
    
    set(water,'P',P,'T',T);
%     set(water_vap,'P',P,'T',T,'X','H2O:1');
    hg = enthalpy_mass(water);
    rhog = density(water);
    ug = intEnergy_mass(water);
    
    setState_Tsat(water,[T 0]);  %liquid
    rhoL = density(water);
    
    %---------------------------------------------
    dT = T - T_prev;
    dQ_metal = m_metal * c_metal * dT;
    dQ_solid = m_solid * c_solid * dT;
    %---------------------------------------------
    dtheta = 0;
    dW = - volume * P * dtheta;
    
    m_a = q_max * m_solid;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P / rho_a;
    du_a = u_a - ua_prev;
    
    dQ_adsorbate = m_a * du_a - dW;
    %----------------------------------------------
    m_g = m_total - m_a;
%     m_g_paper = (volume * density_amb) * (porosity - (q_max * m_solid)/(volume * rho_a1));
    du_g = ug - ug_prev;
    dQ_gas = -dW + m_g * du_g;
    %-----------------------------------------------
    dQ12 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
    %-----------------------------------------------
    Q12 = Q12 + dQ12;
    
    %graphing
    T12_water(i) = T;
    P12_water(i) = P;
    T12_bed(i) = T;
    P12_bed(i) = P;
    q12(i) = q_max;
    m_gas12(i) = m_g;
    m_ads12(i) = m_a;
    Q12_heat(i) = Q12;
    h12(i) = hg;
    m_ref(j) = 0;
    P_ref(j) = P;
    
    ua_prev = u_a;
    T_prev = T;
    i = i + 1;
    j = j + 1;
end


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

Q_heating_ref = hot_m_dot*(hot_h_out - hot_h_in)/1e6  %reference value only


%%--------------------------------------------------------------
%Process 2-3 = Desorption and Condensation;  %constant P
%%--------------------------------------------------------------
% T_cond = cooling_T_in;  %condensation temperature
% setState_Tsat(water,[T_cond 1]); %Vapor phase
% set(water,'T',T_cond,'P',P_cond);
% % P_cond = pressure(water);  %condensation pressure
% hg2 = enthalpy_mass(water);
% ug2 = intEnergy_mass(water);

T_cond = cooling_T_in;  %condensation temperature
setState_Tsat(water,[T_cond 0]); %liquid phase
rhoL = density(water);

setState_Tsat(water,[T_cond,0]);
hf_low = enthalpy_mass(water);

T2_bed = T12_bed(end);
T2_water = T2_bed;
P2_water = P_cond;
P2_bed = P_cond;

T3_water = T_cond;
P3_water = P_cond;
T3_bed = hot_T_in;
P3_bed = P_cond;

% q_min = K*exp(Qst/(R/MolarMass*T_cond))*P_cond;
m_a2 = q_max * m_solid;
m_g2 = m_gas12(end);

%Initial state
set(water,'T',T2_bed,'P',P_cond);
% P_cond = pressure(water);  %condensation pressure
hg2 = enthalpy_mass(water);
ug2 = intEnergy_mass(water);


dT = (T3_bed - T2_bed)/10;
q_prev = q_max;
ha = hg2 - Qst;  %enthalpy of adsorbate phase
rho_a = rhoL;
ua_prev = ha - P_cond / rho_a;
ug_prev = ug2;
ma_prev = m_a2;
mg_prev = m_g2;
T_prev = T2_bed;
Qcond = 0;
Q23 = 0;
Qgas = 0;
m_des = 0;
m_leaving_prev = 0;
m_cond = 0;
i = 1;

for T = T2_bed:dT:T3_bed
    set(water,'T',T,'P',P_cond);
%     set(water_vap,'P',P_cond,'T',T,'X','H2O:1');
    hg = enthalpy_mass(water);
    ug = intEnergy_mass(water);
    rhog = density(water);
    
    mg_capacity = open_volume * rhog * porosity;

    setState_Tsat(water,[T 0]); %liquid
    rhoL = density(water);

    %-----------------------------------------
    dT = T - T_prev;
    dQ_metal = m_metal * c_metal * dT;
    dQ_solid = m_solid * c_solid * dT;
    %-----------------------------------------
    q = K*exp(Qst/(R/MolarMass*T))*P_cond;
    Dq = q - q_prev;
    
    m_a = q * m_solid;
    Dm_a = m_a - ma_prev;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P_cond / rho_a;
    du_a = u_a - ua_prev;
    dU_a = m_a * du_a + u_a * Dm_a;
    
    Dm_des = -Dm_a;
    
    dW_a = - volume * P_cond * Dq;
    
    dQ_adsorbate = dU_a - dW_a - Dm_a * hg;  %changed mass of des to mass of ads
    %---------------------------------------------
%     m_g_paper = (volume * rhog) * (porosity - (q * m_solid)/(volume * rho_a));
    Dm_g = Dm_des;
%     m_g_t = m_total - m_a;

    mg_init = mg_prev + Dm_g; % mass of gas initially
    
    if mg_init > mg_capacity;   % excess mass leaves the adsorber
        m_leaving = mg_init - mg_capacity;
    else
        m_leaving = 0;
    end
    
    m_g = mg_init - m_leaving;
    Dm_g = m_g - mg_prev;
    
    du_g = ug - ug_prev;
    dU_g = m_g * du_g + ug * Dm_g;  %m_g = Dm_a - m_leaving
    
    dW_g = -dW_a;
    
    Dm_cond = m_leaving;
    m_cond = m_cond + Dm_cond;
%     Dm_cond = Dm_g;
%     dq_cond = Dm_cond * (hg);
    dQ_cond = Dm_cond * (hg - hf_low);
    
    
    dQ_gas = dU_g - dW_g - Dm_g * hg; %+ dQ_cond;
    Qgas = Qgas + dQ_gas;
    %---------------------------------------------
    dQ23 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
    %---------------------------------------------

    %---------------------------------------------
    Q23 = Q23 + dQ23;
    Qcond = Qcond + dQ_cond;
    
    m_des = m_des + Dm_des;
    
    %graphing
    T23_water(i) = T;  %?????
    P23_water(i) = P_cond;
    T23_bed(i) = T;
    P23_bed(i) = P_cond;
    q23(i) = q;
    m_desorbed(i) = m_des;
    m_gas23(i) = m_g;
    m_ads23(i) = m_a;
    m_cond23(i) = m_cond;
    Q23_des(i) = Q23;
    Q23_gas(i) = Qgas;
    h23(i) = hg;
    m_ref(j) = m_cond;
    P_ref(j) = P_cond;

    T_prev = T;
    ua_prev = u_a;
    ug_prev = ug;
    q_prev = q;
    ma_prev = m_a;
    mg_prev = m_g;
    m_leaving_prev = m_leaving;
    
    i = i + 1;
    j = j + 1;
end
 
% h_out = cooling_h_in + Q_cond / cooling_m_dot;
% set(water,'P',P_cond,'H',h_out);
T_cooling_out = temperature(water);

Q_cooling_ref = cooling_m_dot*(cooling_h_out - cooling_h_in)/1e6



%State 2
% m_dot_water = cooling_m_dot;
% h2 = h3-Q_cond/m_dot_water;
% P2 = P_cond;
% set(water,'P',P2,'H',h2);
% T2 = temperature(water);

%---------------------------------------------------
%Process 3-4 = Cooling and Depressurization
%----------------------------------------------------
T_max = T3_bed;
P4_bed = P_evap;

% set(water_vap,'P',P_evap,'T',T_amb,'X','H2O:1');
set(water,'T',T_max,'P',P_cond);
ug3 = intEnergy_mass(water);
hg3 = enthalpy_mass(water);
rhog3 = density(water);

setState_Psat(water,[P_evap 0]);  %liquid
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
    dQ_metal = m_metal * c_metal * dT;
    dQ_solid = m_solid * c_solid * dT;
    %---------------------------------------------------
    dtheta = 0; %Va/Vb;
    dW = - volume * P * dtheta;
    
    m_a = q_min * m_solid;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P / rho_a;
    du_a = u_a - ua_prev;
    
%     dQ_adsorbate = -volume * P * dtheta - m_a * du_a;
    dQ_adsorbate = m_a * du_a - dW;
    %---------------------------------------------------
    m_g_paper = (volume * rhog3) * (porosity - (q_min * m_solid)/(volume * rho_a3));
    m_g = m_gas23(end);
    du_g = ug - ug_prev;
    
    dQ_gas = m_g * du_g - dW;
    %---------------------------------------------------
    dQ34 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
    %---------------------------------------------------
    Q34 = Q34 + dQ34;
    
    %graphing
    T34_bed(i) = T;
    P34_bed(i) = P;
    q34(i) = q_min;
    m_gas34(i) = m_g;
    m_ads34(i) = m_a;
    Q34_cool(i) = Q34;
    h34(i) = hg;
    m_ref(j) = m_cond;
    P_ref(j) = P;
    
    T_prev = T;
    ua_prev = u_a;
    ug_prev = ug;
    i = i + 1;
    j = j + 1;
end


%-------------------------------------------------------------
%Process 4-1
%-------------------------------------------------------------

T4_bed = T34_bed(end);
set(water_vap,'P',P_evap,'T',T4_bed,'X','H2O:1');
ug4 = intEnergy_mass(water_vap);
hg4 = enthalpy_mass(water_vap);
rhog4 = density(water_vap);

setState_Tsat(water,[T4_bed 0]); %liquid phase
rhoL = density(water);

m_a4 = q_min * m_solid;
m_g4 = m_gas34(end);

ha = hg4 - Qst;  %enthalpy of adsorbate phase
rho_a = rhoL;
ua_prev = ha - P_cond / rho_a;
q_prev = q_min;
ug_prev = ug4;
ma_prev = m_a4;
mg_prev = m_g4;
T_prev = T4_bed;
Q41 = 0;
m_ads = 0;
m_evap = 0;
Q_evap = 0;
i = 1;
dT = -(T4_bed - T1_bed)/10;

for T = T4_bed:dT:T1_bed
    q = K*exp(Qst/(R/MolarMass*T))*P_evap;
    
    set(water,'T',T,'P',P_evap);
%     set(water_vap,'P',P_evap,'T',T,'X','H2O:1');
    hg = enthalpy_mass(water);
    rhog = density(water);
    ug = intEnergy_mass(water);

    mg_capacity = open_volume * rhog * porosity;
    
    setState_Psat(water,[P_evap,1]);
    h_ev = enthalpy_mass(water);
    
    setState_Tsat(water,[T 0]);
    rhoL = density(water);
    
    %--------------------------------------------------
%     dT = T - T_prev;
    dQ_metal = m_metal * c_metal * dT;
    dQ_solid = m_solid * c_solid * dT;
    %--------------------------------------------------
    Dq = q - q_prev;
    m_a = q * m_solid;
    Dm_a = m_a - ma_prev;
    
    ha = hg - Qst;  %enthalpy of adsorbate phase
    rho_a = rhoL;
    u_a = ha - P_evap / rho_a;
    du_a = u_a - ua_prev;
    dU_a = m_a * du_a + u_a * Dm_a;

    Dm_ads = Dm_a;
    
    dW_a = - volume * P_evap * Dq;
    
    dQ_adsorbate = -h_ev * Dm_ads + dU_a - dW_a;
    %-------------------------------------------------
%     m_g_paper = (volume * density_amb) * (porosity - (q_max * m_solid)/(volume * rho_a1));

   
    m_g = mg_capacity;
    Dm_g = m_g - mg_prev;
    Dm_entering = Dm_ads - Dm_g;
 
    du_g = ug - ug_prev;
    
    dU_g = m_g * du_g + ug * Dm_g;
    
    dQ_gas = dU_g - dW_a - Dm_g * h_ev;
    %--------------------------------------------------
    dQ41 = dQ_metal + dQ_solid + dQ_adsorbate + dQ_gas;
    %--------------------------------------------------
    Dm_evap = Dm_entering;
    m_evap = m_evap + Dm_evap;
    
    dQ_evap = Dm_evap * (h_ev);
    %--------------------------------------------------
    Q41 = Q41 + dQ41;
    Q_evap = Q_evap + dQ_evap;
    
    m_ads = m_ads + Dm_ads;
    
    %graphing
%     T41_water(i) = T;  %?????
    P41_water(i) = P_evap;
    T41_bed(i) = T;
    P41_bed(i) = P_evap;
    q41(i) = q;
    m_adsorbed(i) = m_ads;
    m_gas41(i) = m_g;
    m_ads41(i) = m_a;
    m_evap41(i) = m_evap;
    Q41_ads(i) = Q41;
    h41(i) = hg;
    m_ref(j) = m_cond - m_evap;
    P_ref(j) = P_evap;
  
    T_prev = T;
    q_prev = q;
    ua_prev = u_a;
    ug_prev = ug;
    ma_prev = m_a;
    mg_prev = m_g;
    
    i = i + 1;
    j = j + 1;
end

%Condensation

% % Throttling
% P4_water = P_evap;
% setState_Psat(water,[P_evap,0]); %liquid
% hf = enthalpy_mass(water);
% T_evap = temperature(water);
% T4_water = T_evap;
% 
% setState_Psat(water,[P_evap,1]); %vapor
% hg = enthalpy_mass(water);
% 
% hfg = hf - hg;
% 
% set(water,'T',T_amb,'P',P_evap);
% hf_amb = enthalpy_mass(water);
% hf_evap = hf;
% 
% y = 1 - (hf_amb - hf_evap) / hfg;
% 
% %Evaporation
% Qevap = m_evap * hfg;

%Water states
set(water,'T',T1_water,'P',P1_water);
h1 = enthalpy_mass(water);
set(water,'T',T2_water,'P',P2_water);
h2 = enthalpy_mass(water);
set(water,'T',T3_water,'P',P3_water);
h3 = enthalpy_mass(water);
h4 = h3;

%Throttling
setState_Psat(water,[P_evap,0]); %liquid
hf = enthalpy_mass(water);
setState_Psat(water,[P_evap,1]); %vapor
hg = enthalpy_mass(water);
T4_water = temperature(water);
P4_water = P_evap;
x = (h4 - hf)/(hg - hf);
y = 1-x;

%Evaporation
hfg = hg - hf
Qevap = y * m_evap * (hg - hf);


[h_dome, s_dome, u_dome, T_dome, P_dome, v_dome] = vaporDome(water);

fprintf('\n***************************************************************\n')
m_cond;
m_evap;

Q12 = Q12/1e6
Q23 = Q23/1e6
Q34 = Q34/1e6
Q41 = Q41/1e6
Qevap = Qevap/1e6
Qcond = Qcond/1e6
Qin = (Q12+Q23)

COP = Qevap/Qin

Qadd = (Q12 + Q23 + Qevap)
Qrej = -(Q34 + Q41 - Qcond)


%% Figures

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
% xlabel('-1/T')
% ylabel('lnP')
% legend('bed','water')
% title('Adsorption Chiller Cycle')

% figure(3) % Temperature vs Concentration
% clf
% hold on
% plot(T12_bed-273,q12,'x-')
% plot(T23_bed-273,q23,'x-')
% plot(T34_bed-273,q34,'x-')
% plot(T41_bed-273,q41,'x-')
% xlabel('T (C)')
% ylabel('q (kg/kg)')
% % legend('water','bed')
% title('Adsorption Chiller Concentrations')

% figure(4)
% clf
% hold on
% plot(T12_bed-273,m_ads12,'bx-')
% plot(T23_bed-273,m_ads23,'rx-')
% plot(T34_bed-273,m_ads34,'gx-')
% plot(T41_bed-273,m_ads41,'kx-')
% % plot(T12_bed-273,m_gas12,'bo-')
% % % plot([T12_bed T23_bed]-273,[m_ads12 + m_gas12 m_ads23 + m_gas23],'g-')
% % plot(T23_bed-273,m_gas23,'ro-')
% % plot(T34_bed-273,m_gas34,'go-')
% % plot(T41_bed-273,m_gas41,'ko-')
% xlabel('T (C)')
% ylabel('mass (kg)')
% % legend('adsorbed phase','gas phase','total')
% legend('Heat/Pressure','Desorb/Cond','Cool/Depressure','Evap/Adsorb')
% title('Adsorption Chiller Mass')

figure(5)
clf
hold on
plot(T12_bed-273,m_gas12,'bo-')
plot(T23_bed-273,m_gas23,'ro-')
plot(T34_bed-273,m_gas34,'go-')
plot([T34_bed(end) T41_bed(1)]-273,[m_gas34(end) m_gas41(1)],'ko-')
plot(T41_bed-273,m_gas41,'ko-')
xlabel('T (C)')
ylabel('Vapor Mass (kg)')
legend('Heat/Pressure','Desorb/Cond','Cool/Depressure','Evap/Adsorb')
title('Vapor Mass in Adsorber')



% figure(6)
% clf
% hold on
% plot(T12_bed-273,Q12_heat/1e6,'bx-')
% plot(T23_bed-273,Q23_des/1e6,'bx-')
% xlabel('T (C)')
% ylabel('Driving energy (MJ)')
% title('Adsorption Chiller')
% % axis([45 90 0 180])
% 
% figure(7)
% clf
% hold on
% plot(T34_bed-273,-Q34_cool/1e6,'bx-')
% plot(T41_bed-273,-Q41_ads/1e6,'bx-')
% xlabel('T (C)')
% ylabel('Driving energy (MJ)')
% title('Adsorption Chiller')
% % axis([45 90 0 180])

% figure(8)
% clf
% hold on
% plot(T23_bed-273,Q23_des/1e6,'bx-')
% plot(T23_bed-273,Q23_gas/1e6,'rx-')
% xlabel('T (C)')
% ylabel('Driving energy (MJ)')
% title('Adsorption Chiller')
% legend('Total','Desorbed')
% % axis([45 90 0 180])

figure(9)
clf
hold on
plot([h1 h2 h3 h4 h1]./1e3,...
    [P1_water P2_water P3_water P4_water P1_water]./1e3,'rx-')
plot(h12./1e3,P12_bed./1e3,'bx-')
plot(h23./1e3,P23_bed./1e3,'bx-')
plot(h34./1e3,P34_bed./1e3,'bx-')
plot(h41./1e3,P41_bed./1e3,'bx-')
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

% 
% figure(10)
% clf
% hold on
% plot(T23_bed-273,m_cond23,'bo-')
% plot(T41_bed-273,m_evap41,'go-')
% xlabel('Temperature of Bed (C)')
% ylabel('Refrigerant Mass (kg)')
% legend('Mass Condensed','Mass Evaporated')
% text(T2_bed-273,m_cond23(1),'2')
% text(T3_bed-273,m_cond23(end),'3')
% text(T4_bed-273,m_evap41(1),'4')
% text(T1_bed-273,m_evap41(end),'1')
% title('Refrigerant Mass')

figure(10)
clf
hold on
plot(P_ref/1e3,m_ref,'bo-')
% plot(P23_bed/1e3,m_cond23,'bo-')
% plot(P41_bed/1e3,m_evap41,'go-')
xlabel('Pressure Bed (kPa)')
ylabel('Refrigerant Mass (kg)')
text(P2_bed/1e3,m_cond23(1),'2')
text(P3_bed/1e3,m_cond23(end),'3')
text(P4_bed/1e3,m_evap41(1),'1')
text(P1_bed/1e3,m_evap41(end),'4')
text(3.9,80,'Condensation')
text(3,140,'Throttling')
text(1.7,80,'Evaporation')
text(3,10,'Heating')
title('Refrigerant Mass')

plotfixer