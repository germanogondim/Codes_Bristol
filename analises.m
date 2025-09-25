%% Analises
clc
clear all

% Abertura
load('Extracoes.mat')
data=[data(:,1) data(:,2) ones(size(data,1),1)];

%% 
for iii=1:9
    CC(iii,1)=corr(TEMP_M(:,iii+4),ASW(:,iii+4));
    CC(iii,2)=corr(TEMP_M(:,iii+4),AGW(:,iii+4));
    CC(iii,3)=corr(ATEMP(:,iii+4),AASW(:,iii+4));
    CC(iii,4)=corr(ATEMP(:,iii+4),AAGW(:,iii+4));
end
st=3;
for iii=st:size(ASW,1)
    ST_ASW(iii-st+1,:)=sum(ASW(iii-st+1:iii,:),1);
    ST_AGW(iii-st+1,:)=sum(AGW(iii-st+1:iii,:),1);
    ST_PREC(iii-st+1,:)=sum(PREC_M(iii-st+1:iii,:),1);
end
%% Boxplot Clim
nplotx=3;
nploty=3;
espvbi=0.025;
espvbs=0.035;
espvf=0.085;
esphf=0.020;
esphbe=0.005;
esphbd=0.050;
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
close all
figure('color',[1 1 1],'position',[10 10 1200 600])
c=0;
for iii=1:nplotx
    for ii=1:nploty
        clear agw asw
        c=c+1;
        for i=1:12
            xd=find(data(:,2)==i);
            agw(:,i)=AGW(xd,c+4);
        end
        subplot('position',PLOT{iii,ii})
        boxplot(agw,'symbol','')
        title(num2str(c))
        grid on
    end

end

figure('color',[1 1 1],'position',[10 10 1200 600])
c=0;
for iii=1:nplotx
    for ii=1:nploty
        clear agw asw
        c=c+1;
        for i=1:12
            xd=find(data(:,2)==i);
            asw(:,i)=ASW(xd,c+4);
        end
        subplot('position',PLOT{iii,ii})
        boxplot(asw,'symbol','')
        title(num2str(c+4))
        grid on
    end
end
%% Plot temporal
nplotx=2;
nploty=1;
espvbi=0.065;
espvbs=0.035;
espvf=0.095;
esphf=0.0;
esphbe=0.025;
esphbd=0.050;
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
figure('color',[1 1 1],'position',[10 10 1200 900])

subplot('position',PLOT{1,1})
yyaxis left
plot(datenum(data),TEMP_M(:,6))
yyaxis right
plot(datenum(data),ASW(:,6),'r','LineWidth',1.5),hold on
plot(datenum(data),AGW(:,6),'color',[0.8 0 0.8],'LineWidth',1.5,'LineStyle','-')
grid on
set(gca,'xlim',[datenum(data(1,:)) datenum(data(end,:))],'xtick',[datenum(data(1,:)):180:datenum(data(end,:))])
datetick('x','mm-yyyy','keepticks')
title('02')

subplot('position',PLOT{2,1})
yyaxis left
plot(datenum(data),TEMP_M(:,11))
yyaxis right
plot(datenum(data),ASW(:,11),'r','LineWidth',1.5),hold on
plot(datenum(data),AGW(:,11),'color',[0.8 0 0.8],'LineWidth',1.5,'LineStyle','-')
grid on
set(gca,'xlim',[datenum(data(1,:)) datenum(data(end,:))],'xtick',[datenum(data(1,:)):180:datenum(data(end,:))])
datetick('x','mm-yyyy','keepticks')
title('07')

% subplot('position',PLOT{3,1})
% yyaxis left
% plot(datenum(data),TEMP_M(:,7))
% yyaxis right
% plot(datenum(data),ASW(:,7),'r'),hold on
% plot(datenum(data),AGW(:,7),'color',[0.8 0 0.8])
% grid on
% set(gca,'xlim',[datenum(data(1,:)) datenum(data(end,:))])
% datetick('x','mm-yy','keeplimits')
% title('03')
% 
% subplot('position',PLOT{4,1})
% yyaxis left
% plot(datenum(data),TEMP_M(:,11))
% yyaxis right
% plot(datenum(data),ASW(:,11),'r'),hold on
% plot(datenum(data),AGW(:,11),'color',[0.8 0 0.8])
% grid on
% set(gca,'xlim',[datenum(data(1,:)) datenum(data(end,:))])
% datetick('x','mm-yy','keeplimits')
% title('07')
