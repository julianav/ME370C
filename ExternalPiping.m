function[Hot Chil Cool] = ExternalPiping(hot_T_in,chilled_T_in,cooling_T_in,...
    hot_T_diff,chilled_T_diff,cooling_T_diff)

water = importPhase('liquidVapor.xml','water');
global Po

%Hot water stream
hot_V_dot = 270; %m^3/h
density_water = 996;
hot_m_dot = hot_V_dot*density_water/3600; %kg/s
% hot_T_in = 85 + 273; %K
% hot_T_diff = 5.6; %K
hot_T_out = hot_T_in - hot_T_diff;
hot_P_in = Po; %Pa
hot_P_drop = 0;%49033.25;%29419.95; %Pa = 3mH2O %44100; %Pa = 4.5 mH20
hot_P_out = hot_P_in - hot_P_drop;
set(water,'P',hot_P_in,'T',hot_T_in);
QualityIn = vaporFraction(water);
hot_h_in = enthalpy_mass(water);
HotFlowXin = flowExergy_mass(water);
set(water,'P',hot_P_out,'T',hot_T_out);
QualityOut = vaporFraction(water);
hot_h_out = enthalpy_mass(water);
HotFlowXout = flowExergy_mass(water);

Qheating = hot_m_dot*(hot_h_out - hot_h_in)/1e3; %kW

Hot.Qheating = Qheating;
Hot.h_in = hot_h_in;
Hot.h_out = hot_h_out;
Hot.m_dot = hot_m_dot;
Hot.Flowxin = HotFlowXin;
Hot.Flowxout = HotFlowXout;

%Chilled water stream 
chilled_V_dot = 181; %m^3/h
chilled_m_dot = chilled_V_dot*density_water/3600;
% chilled_T_in = 14 + 273; %K
% chilled_T_diff = 6;
chilled_T_out = chilled_T_in - chilled_T_diff; %K
chilled_P_in = Po; %Pa 
chilled_P_drop = 75500;% Pa  53900; %Pa = 5.5 mH20
chilled_P_out = chilled_P_in - chilled_P_drop;
set(water,'P',chilled_P_in,'T',chilled_T_in);
chilled_h_in = enthalpy_mass(water);
ChilFlowXin = flowExergy_mass(water);
set(water,'P',chilled_P_out,'T',chilled_T_out);
chilled_h_out = enthalpy_mass(water);
ChilFlowXout = flowExergy_mass(water);

Qchilled = chilled_m_dot*(chilled_h_out - chilled_h_in)/1e3;%kW

Chil.Qchilled = Qchilled;
Chil.h_in = chilled_h_in;
Chil.h_out = chilled_h_out;
Chil.m_dot = chilled_m_dot;
Chil.Flowxin = ChilFlowXin;
Chil.Flowxout = ChilFlowXout;


%Cooling water stream
cooling_V_dot = 637; %m^3/h
cooling_m_dot = cooling_V_dot*density_water/3600;
% cooling_T_in = 31 + 273; %K
% cooling_T_diff = 3.8;
cooling_T_out = cooling_T_in + cooling_T_diff;
cooling_P_in = Po; %Pa
cooling_P_drop = 0;%44100; %78453.2; %Pa = 8 mH20
cooling_P_out = cooling_P_in - cooling_P_drop;
set(water,'P',cooling_P_in,'T',cooling_T_in);
cooling_h_in = enthalpy_mass(water);
coolFlowxin = flowExergy_mass(water);
set(water,'P',cooling_P_out,'T',cooling_T_out);
cooling_h_out = enthalpy_mass(water);
coolFlowxout = flowExergy_mass(water);

Qcooling = cooling_m_dot*(cooling_h_out - cooling_h_in)/1e3;%kW

Cool.Qcooling = Qcooling;
Cool.h_in = cooling_h_in;
Cool.h_out = cooling_h_out;
Cool.m_dot = cooling_m_dot;
Cool.Flowxin = coolFlowxin;
Cool.Flowxout = coolFlowxout;

end
