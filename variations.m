clear all
format compact
fprintf('\n***************************************************************\n')
fprintf('NEW RUN')
fprintf('\n***************************************************************\n')

steps = 20;

%Hot water stream
% hot_T_in_range = linspace(50,100,steps)+ 273; %K
 hot_T_in_range = [70 75 80 85 90 95]+273;
 hot_T_in = 65+273;
 hot_T_diff = 5.6; %K

%Chilled water stream 
% chilled_T_out = 9 + 273; %K
chilled_T_in = 12 + 273; %K
chilled_T_diff = 6;

%Cooling water stream
cooling_T_in = 28 + 273; %K
cooling_T_diff = 33-28;

j = 1;
for i = 1%:length(hot_T_in_range) %steps
    fprintf('\n***************************************************************')
    fprintf('\n***************************************************************\n')
%     hot_T_in = hot_T_in_range(i);
%     hot_T_in = 85+273;
    chilled_T_in;
    cooling_T_in;
    
    fprintf('Hot T = %0.f C \n',hot_T_in-273)
    fprintf('Chilled T = %d C \n',chilled_T_in-273)
    fprintf('Cooling T = %d C \n',cooling_T_in-273)
    
    fprintf('***************************************************************\n')
    
    Outputs = multicycle(hot_T_in,chilled_T_in,cooling_T_in,hot_T_diff,chilled_T_diff,cooling_T_diff);
    
    error = Outputs.error;
%     if error < 2 && Outputs.Xadsorb > 0
        COP(j) = Outputs.COP;
        CoolingCap(j) = Outputs.CoolingCap;
        HotT(j) = hot_T_in;
        Xheating(j) = Outputs.Xheating;
        Xdesorb(j) = Outputs.Xdesorb;
        Xcond(j) = Outputs.Xcond;
        Xcooling(j) = Outputs.Xcooling;
        Xevap(j) = Outputs.Xevap;
        Xadsorb(j) = Outputs.Xadsorb;
        Xexpan(j) = Outputs.Xexpan;
        Xloss_system(j) = Outputs.Xsystem;    
        Xefficiency(j) = Outputs.Xefficiency;
        j = j+1;
%     end


end

figure(1)
clf
plot(HotT - 273,COP)
xlabel('Temperature of the Hot Stream (C)')
ylabel('COP')
plotfixer

% figure(2)
% clf
% plot(HotT - 273,CoolingCap)
% xlabel('Temperature of the Hot Stream (C)')
% ylabel('Cooling Capacity kW')
% plotfixer
% 
% figure(3)
% clf
% plot(HotT - 273,Xloss_system)
% xlabel('Temperature of the Hot Stream (C)')
% ylabel('Exergy loss MJ')
% plotfixer
% 
% Ti = 1;
% for i = Ti
%     Xbar = [Xheating(i) Xdesorb(i) Xcond(i) Xcooling(i) Xexpan(i)...
%         Xevap(i) Xadsorb(i)];
% end
% Xlabels = {'1-Heating','2-Desorption','3-Condensation','4-Cooling',...
%     '5-Expansion','6-Evaporation','7-Adsorption'};
% Xloss_total = 0;
% for i=1:1:length(Xbar)
%     Xloss_total = Xloss_total + Xbar(i);
% end
% 
% figure(4)
% clf
% bar(100*Xbar/Xloss_total)
% xlabel('Conversion Device')
% ylabel('Exergy Loss (%)')
% scale = axis;
% text(0.5,56,'Device Key:')
% for i=1:1:length(Xbar)
%     text(0.5,56-4*i,Xlabels(i))
% end
% % text(1,62.5,sprintf('C-C Efficiency:  %3.1f (Exergy),  %3.1f (LHV)',...
% %     100*Exergy_Efficiency,100*Eff_combined_cycle))
% text(0,62.5,sprintf('Hot Temperature: %.2f C  Exergy Efficiency: %.2f',...
%     HotT(Ti)-273,Xefficiency(Ti)))
% axis([0 8 0 60])
% plotfixer
% 
% 
% % for i = 1:length(Xheating)
% %     Xline(i) = [Xheating(i) -Xdesorb(i) Xcond(i) Xcooling(i) Xexpan(i)...
% %         Xevap(i) Xadsorb(i)];
% % end
% Xlabels = {'Heating','Desorption','Condensation','Cooling',...
%     'Expansion','Evaporation','Adsorption'};
% % Xloss_total = 0;
% % for i=1:1:length(Xbar)
% %     Xloss_total = Xloss_total + Xbar(i);
% % end
% 
% figure(5)
% clf
% hold on
% plot(HotT-273,Xheating,'ro-')
% plot(HotT-273,Xdesorb,'bx-')
% plot(HotT-273,Xcond,'gd-')
% plot(HotT-273,Xcooling,'cs-')
% plot(HotT-273,Xexpan,'m+-')
% plot(HotT-273,Xevap,'k*-')
% plot(HotT-273,Xadsorb,'y^-')
% Legend('Heating','Desorption','Condensation','Cooling',...
%     'Expansion','Evaporation','Adsorption')
% xlabel('Temperature of the Hot Stream (C)')
% ylabel('Exergy loss (MJ)')
% plotfixer
% 
% figure(6)
% clf
% hold on
% plot(HotT-273,100*Xheating./Xloss_system,'ro-')
% plot(HotT-273,100*Xdesorb./Xloss_system,'bx-')
% plot(HotT-273,100*Xcond./Xloss_system,'gd-')
% plot(HotT-273,100*Xcooling./Xloss_system,'cs-')
% plot(HotT-273,100*Xexpan./Xloss_system,'m+-')
% plot(HotT-273,100*Xevap./Xloss_system,'k*-')
% plot(HotT-273,100*Xadsorb./Xloss_system,'y^-')
% Legend('Heating','Desorption','Condensation','Cooling',...
%     'Expansion','Evaporation','Adsorption')
% xlabel('Temperature of the Hot Stream (C)')
% ylabel('Exergy loss (%)')
% plotfixer



