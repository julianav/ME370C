function [ T ] = T_isosteric(q,P)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global K Qst open_volume
% Adsorption Model parameters
R = 8314; %J/kmol*K
MolarMass = 18.016;
Rwater = R/MolarMass; %J/kg*K

T_guess = Qst/(Rwater*log(q/(K*P)));

% Ko = 7.3e-10; %kPa-1
% Q = 2693; %kJ/kg 
% Rw = 0.461;
% q_eq = 0.4;
% t = 6;
% y = exp(Q/(Rw*T_guess));
% q_guess = (Ko*y*P/1e3)/(1+((Ko/q_eq)*y*P/1e3)^t)^(1/t);
% q
% 
% V = open_volume;
% n = (P*V)/(R/1e3*T_guess);
% dT = 1;
% dQ = 1;
% % while abs(n - n_init) > .01
%     
%     while abs(q_guess - q) > .01
%         if q_guess > q
%             T_guess = T_guess+dT;
%         elseif q_guess < q
%             T_guess = T_guess-dT;
%         end
%         y = exp(Q/(Rw*T_guess));
%         q_guess = (Ko*y*P/1e3)/(1+((Ko/q_eq)*y*P/1e3)^t)^(1/t);
%     end
% %     n = (P*V)/(R/1e3*T_guess)
% %     Q = Q+dQ;
% %     y = exp(Q/(Rw*T_guess));
% %     q_guess = (Ko*y*P/1e3)/(1+((Ko/q_eq)*y*P/1e3)^t)^(1/t)

T = T_guess;

end

