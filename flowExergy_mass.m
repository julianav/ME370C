function [ x ] = flowExergy_mass(fluid)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

global To Po

h = enthalpy_mass(fluid);
s = entropy_mass(fluid);
b = h - To*s;

ref = importPhase('liquidVapor.xml','water');
set(ref,'T',To,'P',Po);
h_tm = enthalpy_mass(ref);
s_tm = entropy_mass(ref);

g_tm = h_tm - To*s_tm;
% gibbs = gibbs_mass(fluid)

x = b - g_tm;


end

% global m_solid m_metal c_metal c_solid Qst
% ua = Ads.u;
% sa = Ads.s;
% va = Ads.v;
% ma = Ads.m;
% 
% ug = Gas.u;
% sg = Gas.s;
% vg = Gas.v;
% mg = Gas.m;
% 
% uo = dead.uo;
% so = dead.so;
% vo = dead.vo;
% To = dead.To;
% Po = dead.Po;
% 
% DU = (m_metal*c_metal + m_solid*c_solid) * (T-To)...
%     + ma*(ua - uo) + mg*(ug - uo);
% 
% DS = (m_metal*c_metal + m_solid*c_solid) * log(T/To)...
%     + ma*(sa - so) + mg*(sg - uo);
% 
% DV = ma*(va - vo) + mg*(vg - vo);
% 
% A = DU + P*DV - To*DS
% X = A + ma*Qst*(1-To/T)

