function [ Evap ] = Evaporator( hfg,m,xin,xout,Ts )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

global To

Qevap = m * hfg;

Xin = m*xin + Qevap*(1-(To/Ts));
Xout = m*xout;
Xloss = Xin - Xout;


Evap.Qevap = Qevap;
Evap.Xloss = Xloss;


end

% q_max = Adsorbate_Con_Ratio(T_amb,P_evap);
% q_min = Adsorbate_Con_Ratio(T_max,P_cond);
% Qevap = (q_max - q_min) * hfg * m_solid ;
