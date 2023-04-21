function [Tg,CTE,CTE_nonfit,Temperature,TMA] = TMA2Tg_CTE_cursor_fit(DataName,RowNumber,SampleLength,TMA_range,Tg_points)

clc;
figure
%RowNumber is range where data exsists in ACS file.
%FitTempLengt is range where fitting will be conducted.
%The unit of SampleLength must be 'mm'.
%TMA_range = [startTmp endTemp].
%If Tg_points = 0, Tg will be determined from clicked point.

A = importdata(DataName,'\t',RowNumber);

Temperature = A.data(:,2);
TMA = A.data(:,4);

P(1) = plot(Temperature,TMA,'Linewidth',2);

%------------find Tg ------------

if Tg_points == 0
    TemperatureLength = ginput(2);
else
    TemperatureLength = Tg_points;
end

TMALength1 = TMA(find(Temperature > TemperatureLength(1,1)-1 & Temperature < TemperatureLength(1,1)+1));
TemperatureLength1 = Temperature(find(Temperature > TemperatureLength(1,1)-1 & Temperature < TemperatureLength(1,1)+1));

TMALength2 = TMA(find(Temperature > TemperatureLength(2,1)-1 & Temperature < TemperatureLength(2,1)+1));
TemperatureLength2 = Temperature(find(Temperature > TemperatureLength(2,1)-1 & Temperature < TemperatureLength(2,1)+1));

f1 = fit(TemperatureLength1,TMALength1,'poly1');
f2 = fit(TemperatureLength2,TMALength2,'poly1');

hold on
P(2) = plot(f1,':r');
P(3) = plot(f2,':r');
hold off

syms x
Tg = double(solve((f1.p1-f2.p1)*x+(f1.p2-f2.p2) == 0));

hold on
P(4) = plot(Tg,f1.p1*Tg+f1.p2,'+');
plot(TemperatureLength(1,1),f1.p1*TemperatureLength(1,1)+f1.p2,'b.',TemperatureLength(2,1),f2.p1*TemperatureLength(2,1)+f2.p2,'b.')

hold off

xlim([Temperature(1)-10,Temperature(end)+10]);
ylim([min(TMA)-(max(TMA)-min(TMA))/10,max(TMA)+(max(TMA)-min(TMA))/10]);

tg = string(Tg);
text(Tg,f1.p1*Tg+f1.p2-(max(TMA)-min(TMA))/10,tg)

%legend('TMA','fitting','fitting','Tg','Location','Best')
xlabel('Temperature(℃)')
ylabel('TMA/μm')
f = gca;
f.LineWidth = 1.2;

%------------find CTE------------
TMA_fit_Temp = Temperature(find(Temperature > TMA_range(1) & Temperature < TMA_range(2)));
TMA_fit = TMA(find(Temperature > TMA_range(1) & Temperature < TMA_range(2)));

f_tma = fit(TMA_fit_Temp,TMA_fit,'poly1');

hold on
P(5) = plot(f_tma,'b--');
plot(TMA_range,[f_tma.p1*TMA_range(1)+f_tma.p2,f_tma.p1*TMA_range(2)+f_tma.p2],"*b")
hold off 
legend(P([1:5]),{'TMA','fitting','fitting','Tg','CTE fitting'},'Location','Best')

CTE = f_tma.p1/(SampleLength*1000);

CTE_nonfit = (TMA_fit(end) - TMA_fit(1))/(SampleLength*1000)/(TMA_fit_Temp(end) - TMA_fit_Temp(1));

text(sum(TMA_range)/2,(f_tma.p1*TMA_range(1)+f_tma.p2 + f_tma.p1*TMA_range(2)+f_tma.p2)/2,string(CTE));
xlabel("Temperature (℃)",'FontName','Times','FontSize',15)
ylabel("Displacement /\mu m",'FontName','Times','FontSize',15)
