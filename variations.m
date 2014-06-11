clear all
close all
format compact
fprintf('\n***************************************************************\n')
fprintf('NEW RUN')
fprintf('\n***************************************************************\n')

%Efficiency Modes
% 0 = default, with no efficiciency mode
% 1 = heat recovery mode for steady state conditions
EffMode = 0;

TempMode = 2;
switch TempMode
    case 1 %Standard operating conditions
        hot_T_in_range = 85+273;
        cooling_T_in_range = 31 + 273; %K  
    case 2 %Hot temp range, cool SOC
        hot_T_in_range = [55 60 65 70 75 80 85 90 95 100]+273;
        cooling_T_in_range = 31 + 273; %K 
    case 3 %Hot and cool temp range
        hot_T_in_range = [55 60 65 70 75 80 85 90 95 100]+273;
        cooling_T_in_range = [25 30 35]+273;
end

%Hot water stream
hot_T_diff = 5.6; %K

%Chilled water stream
chilled_T_in = 14 + 273; %K %Standard operating conditions
chilled_T_diff = 5;

%Cooling water stream
cooling_T_diff = 3.8;


for m = 1:length(cooling_T_in_range)
    cooling_T_in = cooling_T_in_range(m);    

    k = 1;
for i = 1:length(hot_T_in_range)
    fprintf('\n***************************************************************')
    fprintf('\n***************************************************************\n')
    hot_T_in = hot_T_in_range(i);
    chilled_T_in;
      
    fprintf('Hot T = %0.f C \n',hot_T_in-273)
    fprintf('Chilled T = %0.f C \n',chilled_T_in-273)
    fprintf('Cooling T = %0.1f C \n',cooling_T_in-273)
    
    fprintf('***************************************************************\n')
    
    Outputs = CrossOver(TempMode,EffMode,hot_T_in,chilled_T_in,cooling_T_in,hot_T_diff,chilled_T_diff,cooling_T_diff);
    
    if Outputs.Qdotchil >0
        %Outputs.Xviolation ~= 1  %checks to make sure the unit is in operating range
        TempLift(m,k) = hot_T_in - cooling_T_in; %Generation temperature lift
        COP(m,k) = Outputs.COP;
        CoolingCap(m,k) = Outputs.CoolingCap;
        HotT(m,k) = hot_T_in;
        CoolT(m,k) = cooling_T_in;
        Xheating(m,k) = Outputs.Xheating;
        Xcooling(m,k) = Outputs.Xcooling;
        Xchil(m,k) = Outputs.Xchill;
        Xexpan(m,k) = Outputs.Xexpan;
        Xsystem(m,k) = Outputs.Xsystem;
        Xefficiency(m,k) = Outputs.Xefficiency;
        Qheating(m,k) = Outputs.Qdothot;
        Qevap(m,k) = Outputs.Qdotchil;
        qmax = Outputs.qmax;
        qmin = Outputs.qmin;
        Pheating = Outputs.Pheating;
        m_ref(m,k) = Outputs.mr;
        
    else
       TempLift(m,k) = nan;
       COP(m,k) = nan;
       CoolingCap(m,k) = nan;
       HotT(m,k) = nan;
    end

    k = k+1;
end

end

%-----------------------------------------------------------
% FIGURES
%-----------------------------------------------------------

if TempMode == 1 %SOC
    Ti = 1;
    for i = Ti
        Xbar = [Xheating(i) Xcooling(i) Xexpan(i) Xchil(i)];
    end
    Xlabels = {'1-Heating/Desorption','2-Cooling/Adsorption','3-Expansion','4-Evaporation'};
    Xloss_total = 0;
    for i=1:1:length(Xbar)
        Xloss_total = Xloss_total + Xbar(i);
    end
    figure()
    clf
    bar(100*Xbar./Xloss_total)
    xlabel('Conversion Device')
    ylabel('Exergy Destruction (% of Total Dest)')
    scale = axis;
    text(2.5,66,'Device Key:')
    for i=1:1:length(Xbar)
        text(2.5,66-5*i,Xlabels(i))
    end
    text(0,82.5,sprintf('Hot Temperature: %.0f C  Exergy Efficiency: %.2f',...
        HotT(Ti)-273,Xefficiency(Ti)))
%     axis([0 5 0 80])
    plotfixer
    

elseif TempMode == 2 %only hot range
    figure()
    clf
    plot(HotT - 273,COP,'x-')
    xlabel('Temperature of the Hot Stream (\circC)')
    ylabel('COP')
    axis([50 100 0.5 .7])
    plotfixer
    
    figure()
    clf
    plot(HotT' - 273,CoolingCap','x-')
    xlabel('Temperature of the Hot Stream (\circC)')
    ylabel('Cooling Capacity (kW)')
    plotfixer
    
    figure()
    clf
    plot(HotT - 273,Xsystem/1e3,'x-')
    xlabel('Temperature of the Hot Stream (\circC)')
    ylabel('Exergy Destruction (kW)')
    plotfixer
    
    figure()
    clf
    hold on
    plot(HotT-273,Xheating/1e3,'ro-')
    plot(HotT-273,Xcooling/1e3,'cs-')
    plot(HotT-273,Xchil/1e3,'g*-')
    plot(HotT-273,Xexpan/1e3,'bd-')
    % plot([HotT(1) HotT(end)]-273,[0 0],'k-')
    legend('Heating/Desorption','Cooling/Adsorption',...
        'Evaporation','Expansion')
    xlabel('Temperature of the Hot Stream (\circC)')
    ylabel('Exergy Destruction (kW)')
    plotfixer
    
    figure()
    clf
    hold on
    plot(HotT-273,100*Xheating./Xsystem,'ro-')
    plot(HotT-273,100*Xcooling./Xsystem,'cs-')
    plot(HotT-273,100*Xchil./Xsystem,'g*-')
    plot(HotT-273,100*Xexpan./Xsystem,'bd-')
    legend('Heating/Desorption','Cooling/Adsorption',...
        'Evaporation','Expansion')
    xlabel('Temperature of the Hot Stream (\circC)')
    ylabel('Exergy loss (%)')
    plotfixer
    
    figure()
    hold on
    plot(HotT-273,Xefficiency,'o-')
    xlabel('Temperature of the Hot Stream (\circC)')
    ylabel('Exergy Efficiency')
    plotfixer
    
elseif TempMode == 3 %both temp ranges
    figure()
    clf
    hold on
    plot(TempLift',CoolingCap','x-')
    xlabel('Generation Temperature Lift (HotT - CoolT)')
    ylabel('Cooling Capacity (kW)')
    labels = {'25\circC' '30\circC' '35\circC' '40\circC'};
    legend(labels)
    Tticks = hot_T_in_range(1)-cooling_T_in_range(end):5:hot_T_in_range(end)-cooling_T_in_range(1);
    set(gca,'XTick',Tticks)
    set(gca,'XTickLabel',{Tticks})
    plotfixer
    
    figure()
    clf
    hold on
    plot(transpose(TempLift),transpose(COP),'x-')
    xlabel('Generation Temperature Lift (HotT - CoolT)')
    ylabel('COP')
    coolTlabels = num2str(cooling_T_in_range - 273);
    labels = {'25\circC' '30\circC' '35\circC' '40\circC'};
    legend(labels)
    Tticks = hot_T_in_range(1)-cooling_T_in_range(end):5:hot_T_in_range(end)-cooling_T_in_range(1);
    set(gca,'XTick',Tticks)
    set(gca,'XTickLabel',{Tticks})
    plotfixer

    figure()
    clf
    hold on
    plot(HotT'-273,CoolingCap','x-')
    xlabel('Hot Temperature Inlet (\circC)')
    ylabel('Cooling Capacity (kW)')
    labels = {'25\circC' '30\circC' '35\circC'};
    legend(labels)
    Tticks = hot_T_in_range-273;
    set(gca,'XTick',Tticks)
    set(gca,'XTickLabel',{Tticks})
    plotfixer
    
    figure()
    clf
    hold on
    plot(HotT'-273,COP','x-')
    xlabel('Hot Temperature Inlet (\circC)')
    ylabel('COP')
    coolTlabels = num2str(cooling_T_in_range - 273);
    labels = {'25\circC' '30\circC' '35\circC'};
    legend(labels)
    Tticks = hot_T_in_range-273;
    set(gca,'XTick',Tticks)
    set(gca,'XTickLabel',{Tticks})
    plotfixer

end







