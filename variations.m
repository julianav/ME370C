clear all
format compact
fprintf('\n***************************************************************\n')
fprintf('NEW RUN')
fprintf('\n***************************************************************\n')

steps = 20;

%Hot water stream
hot_T_in_range = linspace(50,100,steps)+ 273; %K

%Chilled water stream 
chilled_T_out = 9 + 273; %K

%Cooling water stream
cooling_T_in = 31 + 273; %K

j = 1;
for i = 1:steps
    fprintf('\n***************************************************************\n')
    fprintf('\n***************************************************************\n')
    hot_T_in = hot_T_in_range(i)%85+273
    chilled_T_out
    cooling_T_in
    
    Outputs = multicycle(hot_T_in,chilled_T_out,cooling_T_in);
    
    error = Outputs.error;
    if error < 2
        COP(j) = Outputs.COP;
        CoolingCap(j) = Outputs.CoolingCap;
        HotT(j) = hot_T_in;
        Xheating(j) = Outputs.Xheating/1e6;
        Xdesorb(j) = Outputs.Xdesorb/1e6;
        Xcond(j) = Outputs.Xcond/1e6;
        Xcooling(j) = Outputs.Xcooling/1e6;
        Xevap(j) = Outputs.Xevap/1e6;
        Xadsorb(j) = Outputs.Xadsorb/1e6;
        Xexpan(j) = Outputs.Xexpan/1e6;
        Xloss_system(j) = Outputs.Xsystem/1e6;     
        j = j+1;
    end


end

figure(1)
clf
plot(HotT - 273,COP)
xlabel('Temperature of the Hot Stream (C)')
ylabel('COP')
plotfixer

figure(2)
clf
plot(HotT - 273,CoolingCap)
xlabel('Temperature of the Hot Stream (C)')
ylabel('Cooling Capacity kW')
plotfixer

figure(3)
clf
plot(HotT - 273,Xloss_system)
xlabel('Temperature of the Hot Stream (C)')
ylabel('Exergy loss MJ')
plotfixer

for i = 13
    Xbar = [Xheating(i) -Xdesorb(i) Xcond(i) Xcooling(i) Xexpan(i)...
        Xevap(i) Xadsorb(i)];
end
Xlabels = {'1-Heating','2-Desorption','3-Condensation','4-Cooling',...
    '5-Expansion','6-Evaporation','7-Adsorption'};
Xloss_total = 0;
for i=1:1:length(Xbar)
    Xloss_total = Xloss_total + Xbar(i);
end

figure(4)
clf
bar(100*Xbar/Xloss_total)
xlabel('Conversion Device')
ylabel('Exergy Loss (%)')
scale = axis;
text(0.5,56,'Device Key:')
for i=1:1:length(Xbar)
    text(0.5,56-4*i,Xlabels(i))
end
% text(1,62.5,sprintf('C-C Efficiency:  %3.1f (Exergy),  %3.1f (LHV)',...
%     100*Exergy_Efficiency,100*Eff_combined_cycle))
axis([0 8 0 60])
plotfixer


% for i = 1:length(Xheating)
%     Xline(i) = [Xheating(i) -Xdesorb(i) Xcond(i) Xcooling(i) Xexpan(i)...
%         Xevap(i) Xadsorb(i)];
% end
Xlabels = {'Heating','Desorption','Condensation','Cooling',...
    'Expansion','Evaporation','Adsorption'};
% Xloss_total = 0;
% for i=1:1:length(Xbar)
%     Xloss_total = Xloss_total + Xbar(i);
% end

figure(5)
clf
hold on
plot(HotT-273,100*Xheating./Xloss_total,'ro-')
plot(HotT-273,100*Xdesorb./Xloss_total,'bx-')
plot(HotT-273,100*Xcond./Xloss_total,'gd-')
plot(HotT-273,100*Xcooling./Xloss_total,'cs-')
plot(HotT-273,100*Xexpan./Xloss_total,'m+-')
plot(HotT-273,100*Xevap./Xloss_total,'k*-')
plot(HotT-273,100*Xadsorb./Xloss_total,'y^-')
Legend('Heating','Desorption','Condensation','Cooling',...
    'Expansion','Evaporation','Adsorption')


