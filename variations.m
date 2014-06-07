clear all
format compact
fprintf('\n***************************************************************\n')
fprintf('NEW RUN')
fprintf('\n***************************************************************\n')

steps = 20;

global P_cond P_evap

%Hot water stream
% hot_T_in_range = linspace(50,100,steps)+ 273; %K
%  hot_T_in_range = [50 55 65 70 75 80 85 90 95 100]+273;
 hot_T_in_range = 88.5+273;
 hot_T_diff = 5.6; %K

%Chilled water stream 
% chilled_T_out = 9 + 273; %K
chilled_T_in = 11 + 273; %K
chilled_T_diff = 5;

%Cooling water stream
cooling_T_in_range = 29.5 + 273; %K
% cooling_T_in_range = [28 29 30 31 32]+273;
cooling_T_diff = 3.8;

j = 1;
for i = 1:length(hot_T_in_range) %steps
    fprintf('\n***************************************************************')
    fprintf('\n***************************************************************\n')
    hot_T_in = hot_T_in_range(i);
    chilled_T_in;
    cooling_T_in = cooling_T_in_range;
    
    fprintf('Hot T = %0.f C \n',hot_T_in-273)
    fprintf('Chilled T = %0.f C \n',chilled_T_in-273)
    fprintf('Cooling T = %0.1f C \n',cooling_T_in-273)
    
    fprintf('***************************************************************\n')
    
    Outputs = multicycle(i,hot_T_in,chilled_T_in,cooling_T_in,hot_T_diff,chilled_T_diff,cooling_T_diff);
    
    error = Outputs.error;
%     if Outputs.Xcooling > 0
        COP(j) = Outputs.COP;
        CoolingCap(j) = Outputs.CoolingCap;
        HotT(j) = hot_T_in;
        CoolT(j) = cooling_T_in;
        Xheating(j) = Outputs.Xheating;
        Xcooling(j) = Outputs.Xcooling;
        Xchil(j) = Outputs.Xchill;
        Xexpan(j) = Outputs.Xexpan;
        Xsystem(j) = Outputs.Xsystem;
        Xefficiency(j) = Outputs.Xefficiency;
        Qheating(j) = Outputs.Qhot;
        Qevap(j) = Outputs.Qchil;
        qmax = Outputs.qmax;
        qmin = Outputs.qmin;
        Pheating = Outputs.Pheating;
        m_ref(j) = Outputs.mr;
        j = j+1;
%     end


end

CInc = (Qevap(end)-Qevap(1))/1e3;
HInc = (Qheating(end)-Qheating(1))/1e3;

% figure(1)
% clf
% plot(CoolT - 273,COP)
% xlabel('Temperature of the Cool Stream (C)')
% ylabel('COP')
% plotfixer
% 
% 
% figure(1)
% clf
% plot(HotT - 273,COP)
% xlabel('Temperature of the Hot Stream (C)')
% ylabel('COP')
% plotfixer
% 
% figure(2)
% clf
% plot(HotT - 273,CoolingCap)
% xlabel('Temperature of the Hot Stream (C)')
% ylabel('Cooling Capacity (kW)')
% plotfixer
% 
% figure(3)
% clf
% plot(HotT - 273,Xsystem/1e3)
% xlabel('Temperature of the Hot Stream (C)')
% ylabel('Exergy Distruction (kW)')
% plotfixer
% 
% Ti = 3;
% for i = Ti
%     Xbar = [Xheating(i) Xcooling(i) Xexpan(i) Xchil(i)];
% end
% Xlabels = {'1-Heating/Desorption','2-Cooling/Adsorption','3-Expansion','4-Evaporation'};
% Xloss_total = 0;
% for i=1:1:length(Xbar)
%     Xloss_total = Xloss_total + Xbar(i);
% end
% 
% figure(4)
% clf
% bar(100*Xbar./Xloss_total)
% xlabel('Conversion Device')
% ylabel('Exergy Destruction (%)')
% scale = axis;
% text(2.5,66,'Device Key:')
% for i=1:1:length(Xbar)
%     text(2.5,66-5*i,Xlabels(i))
% end
% % text(1,62.5,sprintf('C-C Efficiency:  %3.1f (Exergy),  %3.1f (LHV)',...
% %     100*Exergy_Efficiency,100*Eff_combined_cycle))
% text(0,82.5,sprintf('Hot Temperature: %.0f C  Exergy Efficiency: %.2f',...
%     HotT(Ti)-273,Xefficiency(Ti)))
% axis([0 5 0 80])
% plotfixer
% 
% figure(5)
% clf
% hold on
% plot(HotT-273,Xheating/1e3,'ro-')
% plot(HotT-273,Xcooling/1e3,'cs-')
% plot(HotT-273,Xchil/1e3,'g*-')
% plot(HotT-273,Xexpan/1e3,'bd-')
% % plot([HotT(1) HotT(end)]-273,[0 0],'k-')
% legend('Heating/Desorption','Cooling/Adsorption',...
%     'Evaporation','Expansion')
% xlabel('Temperature of the Hot Stream (C)')
% ylabel('Exergy Destruction (kW)')
% plotfixer
% 
% figure(6)
% clf
% hold on
% plot(HotT-273,100*Xheating./Xsystem,'ro-')
% plot(HotT-273,100*Xcooling./Xsystem,'cs-')
% plot(HotT-273,100*Xchil./Xsystem,'g*-')
% plot(HotT-273,100*Xexpan./Xsystem,'bd-')
% legend('Heating/Desorption','Cooling/Adsorption',...
%     'Evaporation','Expansion')
% xlabel('Temperature of the Hot Stream (C)')
% ylabel('Exergy loss (%)')
% plotfixer



