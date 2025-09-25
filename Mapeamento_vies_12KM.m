clc
clear all
%% Abertura
load('pol_12km')
load('UKCP18_12KM.mat')
load('HadUK_Grid_12KM.mat')

% Dates month
c=0;
for ano=1981:2010
    for mes=1:12
        c=c+1;
        datas_m(c,:)=[ano mes];
    end
end
%% Climatology
for mes=1:12
    xd=find(datas_m(:,2)==mes);
    CLIM_OBS12K(mes,:)=mean(PM_OBS12KM(xd,:),1);
end
%% Bias correction factor
for iii=1:12
    xd=find(datas_m(:,2)==iii);
    for ii=1:size(CLIM_S12KM,1)
        FC_12KM{ii,1}(iii,:)=CLIM_OBS12K(iii,:)./CLIM_S12KM{ii,1}(iii,:);
    end
end
%% Smothing

for iii=1:size(FC_12KM,1)
    for ii=1:size(FC_12KM{1,1},1)
        MB=zeros(size(LAT12K))/0;
        MBS=zeros(size(LAT12K))/0;
        MB(in_12k)=FC_12KM{iii,1}(ii,:);
        [x,y]=find(MB>=0);
        for i=1:length(x)
            xx=[x(i)-1:x(i)+1];xx=xx'*ones(1,length(xx)); % 3X3 matrix moving box
            yy=[y(i)-1:y(i)+1];yy=ones(length(yy),1)*yy;  % 3X3 matrix moving box
            xx(xx<=0)=NaN;xx(xx>size(LAT12K,1))=NaN;         % NaN for values outside the boundaries of the Grid
            yy(yy>size(LAT12K,2))=NaN;yy(yy<=0)=NaN;         % NaN for values outside the boundaries of the Grid
            xx=reshape(xx,[9 1]);
            yy=reshape(yy,[9 1]);
            PP=((size(LAT12K,1))*(yy-1))+xx;                 % Moving box position
            mb=zeros(size(PP))/0;
            mb((isnan(PP)==0))=MB(PP(isnan(PP)==0));      % Select only the positions that are inside the boundaries of the Grid
            ww=zeros(size(PP))/0;
            ww(isnan(mb)==0)=0.5/(sum(isnan(mb)==0)-1);   % Weight division according to Guillod et al. (2018).
            ww(5)=0.5;                                    % The center has 0.5 Weight according to Guillod et al. (2018).
            mb=mb.*ww;
            MBS(x(i),y(i))=sum(mb(mb>0));
        end
        ms=MBS(in_12k);
        FCS_12KM{iii,1}(ii,:)=ms;
    end
end
%% PLOT
ff=[1 4 5 6 7 8 9 10 11 12 13 15 23 25 27 29];
for iii=1:length(ff)
    FF{iii,1}=num2str(ff(iii));
end
nplotx=12;
nploty=16;
CC=1:(nplotx*nploty);
CC=reshape(CC,[nplotx nploty]);
espvbi=0.045;
espvbs=0.03;
espvf=0.015;
esphf=0.00;
esphbe=0.025;
esphbd=0.015;
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
figure('color',[1 1 1 1],'position',[100 100 1600 1200])
MB=zeros(size(LAT12K));
MB(in_12k)=1;
[x,y]=find(MB==1);
xmin=min(x);xmax=max(x);
ymin=min(y);ymax=max(y);
for iii=1:nplotx
    for ii=1:nploty
        subplot('position',PLOT{iii,ii})
        MB=zeros(size(LAT12K))/0;
        MB(in_12k)=FCS_12KM{ii,1}(iii,:);
        imagesc(MB(xmin:xmax,ymin:ymax),[0.4 1.6])
        colormap([[1 1 1]; flipud(interp_cmap)])
        axis off
    end
end

f = gcf;
exportgraphics(f,['FC_12KM.png'],'Resolution',300)
%% BOXPLOT
nplotx=4;
nploty=3;
CC=1:(nplotx*nploty);
CC=reshape(CC,[nplotx nploty]);
espvbi=0.035;
espvbs=0.03;
espvf=0.055;
esphf=0.015;
esphbe=0.005;
esphbd=0.045;
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
figure('color',[1 1 1 1],'position',[100 100 1500 700])
for ii=1:12
    for iii=1:16
    BB(:,iii)=FCS_12KM{iii,1}(ii,:);
    end
    [x,y]=find(CC==ii);
    subplot('position',PLOT{x,y})
    boxplot(BB,FF,'symbol','')
    ylim([0.3 2.8])
    title(num2str(ii))
    grid on
end
f = gcf;
exportgraphics(f,['Boxplot_FC_12KM.png'],'Resolution',300)