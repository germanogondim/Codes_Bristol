%% Analise Abstractions
clc
clear all

%% Abertura
load('Abstractions.mat')
load('Extracoes.mat','data')
%% Grid
load('BNG_convert.mat')
LAT_BNG=flipud(double(x_bng)*ones(1,length(y_bng)));
LON_BNG=ones(size(x_bng))*double(y_bng)';
LAT=rot90(ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\Water Abstractions\Data\Gridded_Actual_Abstraction_Data_For_UoB\Gridded_Actual_Abstraction_Data_For_UoB\England_Monthly_Weighted_Abstractions_1km_Grid_SW_199901_201412.nc','lat'));
LON=rot90(ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\Water Abstractions\Data\Gridded_Actual_Abstraction_Data_For_UoB\Gridded_Actual_Abstraction_Data_For_UoB\England_Monthly_Weighted_Abstractions_1km_Grid_SW_199901_201412.nc','lon'));
nlat=size(LAT,1);nlon=size(LON,2);
reg_uk=shaperead('C:\Users\kl23221\OneDrive - University of Bristol\Data\Shapes\regions_outlines_WGS84_london_merged\regions_outlines_WGS84_london_merged.shp');
MB=zeros(size(LAT));
for iii=1:size(reg_uk,1)
    in=inpolygon(LAT(1:nlat*nlon),LON(1:nlat*nlon),reg_uk(iii,1).Y,reg_uk(iii,1).X); % Check points inside polyon
    in=find(in);
    MB(in)=iii;
    IN{iii,1}=in;
    IN{iii,2}=reg_uk(iii,1).geo_region;
end
%% Camels
load('C:\Users\kl23221\OneDrive - University of Bristol\Data\CAMELS-UK\8344e4f3-d2ea-44f5-8afa-86d2987543a9\data\timeseries\topgraphic_attributes')
arq=dir('C:\Users\kl23221\OneDrive - University of Bristol\Data\CAMELS-UK\8344e4f3-d2ea-44f5-8afa-86d2987543a9\data\timeseries\*.csv');
addpath('C:\Users\kl23221\OneDrive - University of Bristol\Data\CAMELS-UK\8344e4f3-d2ea-44f5-8afa-86d2987543a9\data\timeseries\')  
datas_camel=datenum([1970 10 01]):datenum([2015 09 30]);
datas_camel=datevec(datas_camel);
xi=find(datas_camel(:,1)==1999 & datas_camel(:,2)==1 & datas_camel(:,3)==1);
xf=find(datas_camel(:,1)==2014 & datas_camel(:,2)==12 & datas_camel(:,3)==31);
datas_camel=datas_camel(xi:xf,:);
for k=1:size(arq,1)
    camel=camel_read(['C:\Users\kl23221\OneDrive - University of Bristol\Data\CAMELS-UK\8344e4f3-d2ea-44f5-8afa-86d2987543a9\data\timeseries\' arq(k,1).name]);
    temp(:,k)=camel(xi:xf,3);
    prec(:,k)=camel(xi:xf,1);
    vazao(:,k)=camel(xi:xf,5);
    lat=topgraphic_attributes(k,3);
    lon=topgraphic_attributes(k,4);
    [x,~]=find(abs(LAT-lat)==min(min(abs(LAT-lat))));
    [~,y]=find(abs(LON-lon)==min(min(abs(LON-lon))));
    regreg(k,1)=MB(x(1),y(1));
end
%% Analise por regiao
for iii=1:size(reg_uk,1)
    xd=find(regreg==iii);
    tt=temp(:,xd);
    nn=sum(isnan(tt),2);
    tt(isnan(tt))=0;
    tt=sum(tt,2)./(size(tt,2)-nn);
    TEMP(:,iii)=tt;

    tt=prec(:,xd);
    nn=sum(isnan(tt),2);
    tt(isnan(tt))=0;
    tt=sum(tt,2)./(size(tt,2)-nn);
    PREC(:,iii)=tt;

    tt=vazao(:,xd);
    nn=sum(isnan(tt),2);
    tt(isnan(tt))=0;
    tt=sum(tt,2)./(size(tt,2)-nn);
    VAZAO(:,iii)=tt;
end
c=0;
for ano=1999:2014
    for mes=1:12
        xd=find(datas_camel(:,1)==ano & datas_camel(:,2)==mes);
        c=c+1;
        temp_m(c,:)=mean(TEMP(xd,:),1);
        prec_m(c,:)=mean(PREC(xd,:),1);
        vazao_m(c,:)=mean(VAZAO(xd,:),1);
    end
end
%% Separaçao das Abstrações
lfd_medium=[10 20 30 40 50 110 140 160 190 210 260 270 330 340 350 360 380 390 470 ];
lfd_high=[60 80 120 150 280 400 410 420 600 610 620 670];
lfd_verylow=[70 90 100 180 200 220 230 240 250 290 320 370 430 440 450 480 630 640 650 660];
lfd_low=[130 170 300 310 460];
lfd_nonchar=490;

%% Analise por regiao por ano
reg=unique(unique(MB));
reg=reg(reg>0);
REG_GW=zeros(length(1999:2014),6,length(reg));
REG_SW=zeros(length(1999:2014),6,length(reg));
REG_TW=zeros(length(1999:2014),6,length(reg));

REG_GW_M=zeros(size(data,1),6,length(reg));
REG_SW_M=zeros(size(data,1),6,length(reg));
REG_TW_M=zeros(size(data,1),6,length(reg));
%% Abstraçoes mensais
for iii=1:size(XY,1)
    [x,~]=find(abs(LAT_BNG-XY(iii,2))==min(min(abs(LAT_BNG-XY(iii,2)))));x=x(1);
    [~,y]=find(abs(LON_BNG-XY(iii,1))==min(min(abs(LON_BNG-XY(iii,1)))));y=y(1);
    % XY_GRID(iii,:)=[x y];
    mb=MB(x,y);
    mb=find(reg==mb);
    if strcmp(Abstraction_ID{iii,1},'A')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
            REG_GW_M(:,1,mb)=REG_GW_M(:,1,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'M')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
          REG_GW_M(:,2,mb)=REG_GW_M(:,2,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'E')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
        REG_GW_M(:,3,mb)=REG_GW_M(:,3,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'I')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
         REG_GW_M(:,4,mb)=REG_GW_M(:,4,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'P')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
          REG_GW_M(:,5,mb)=REG_GW_M(:,5,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'W')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
          REG_GW_M(:,6,mb)=REG_GW_M(:,6,mb)+(Abstraction(iii,:))';
    elseif  strcmp(Abstraction_ID{iii,1},'A')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
          REG_SW_M(:,1,mb)=REG_SW_M(:,1,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'M')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
          REG_SW_M(:,2,mb)=REG_SW_M(:,2,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'E')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
          REG_SW_M(:,3,mb)=REG_SW_M(:,3,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'I')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
          REG_SW_M(:,4,mb)=REG_SW_M(:,4,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'P')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
          REG_SW_M(:,5,mb)=REG_SW_M(:,5,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'W')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
          REG_SW_M(:,6,mb)=REG_SW_M(:,6,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'A')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
          REG_TW_M(:,1,mb)=REG_TW_M(:,1,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'M')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
          REG_TW_M(:,2,mb)=REG_TW_M(:,2,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'E')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
          REG_TW_M(:,3,mb)=REG_TW_M(:,3,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'I')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
          REG_TW_M(:,4,mb)=REG_TW_M(:,4,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'P')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
          REG_TW_M(:,5,mb)=REG_TW_M(:,5,mb)+(Abstraction(iii,:))';
    elseif strcmp(Abstraction_ID{iii,1},'W')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
          REG_TW_M(:,6,mb)=REG_TW_M(:,6,mb)+(Abstraction(iii,:))';
    end
end


%% Analise de correlação
for iii=4:size(REG_GW_M,3)
    for ii=1:size(REG_GW_M,2)

        corr_temp(iii-3,1,ii)=corr(temp_m(:,iii),REG_SW_M(:,ii,iii),'Type','Spearman','Rows','complete');
        corr_temp(iii-3,2,ii)=corr(temp_m(:,iii),REG_GW_M(:,ii,iii),'Type','Spearman','Rows','complete');
        corr_temp(iii-3,3,ii)=corr(temp_m(:,iii),REG_TW_M(:,ii,iii),'Type','Spearman','Rows','complete');

        corr_prec(iii-3,1,ii)=corr(prec_m(:,iii),REG_SW_M(:,ii,iii),'Type','Spearman','Rows','complete');
        corr_prec(iii-3,2,ii)=corr(prec_m(:,iii),REG_GW_M(:,ii,iii),'Type','Spearman','Rows','complete');
        corr_prec(iii-3,3,ii)=corr(prec_m(:,iii),REG_TW_M(:,ii,iii),'Type','Spearman','Rows','complete');

        corr_vazao(iii-3,1,ii)=corr(vazao_m(:,iii),REG_SW_M(:,ii,iii),'Type','Spearman','Rows','complete');
        corr_vazao(iii-3,2,ii)=corr(vazao_m(:,iii),REG_GW_M(:,ii,iii),'Type','Spearman','Rows','complete');
        corr_vazao(iii-3,3,ii)=corr(vazao_m(:,iii),REG_TW_M(:,ii,iii),'Type','Spearman','Rows','complete');
    end
end
%% PLOT barras
nplotx=1;
nploty=3;
espvbi=0.055;
espvbs=0.045;
espvf=0.035;
esphf=0.050;
esphbe=0.05;
esphbd=0.080;
rg=[0 8];
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
c=0;
TT{1,1}='Agriculture';
TT{2,1}='Amenity';
TT{3,1}='Environmental';
TT{4,1}='Industrial, Commercial and Public Service';
TT{5,1}='Production of Energy';
TT{6,1}='Water Supply';

for iii=1:size(corr_temp,3)
    figure('color',[1 1 1],'position',[10 10 1300 400])
    subplot('position',PLOT{1,1})
    bar([5:13],corr_temp(:,:,iii))
    ylim([-0.8 0.8])
    grid on
    subplot('position',PLOT{1,2})
    bar([5:13],corr_prec(:,:,iii))
    ylim([-0.8 0.8])
    grid on
    subplot('position',PLOT{1,3})
    bar([5:13],corr_vazao(:,:,iii))
    ylim([-0.8 0.8])
    grid on
    f = gcf;
    exportgraphics(f,['corr_' TT{iii,1} '.png'],'Resolution',300)
end
%% Abstraçoes anuais
for iii=1:size(XY,1)
    [x,~]=find(abs(LAT_BNG-XY(iii,2))==min(min(abs(LAT_BNG-XY(iii,2)))));x=x(1);
    [~,y]=find(abs(LON_BNG-XY(iii,1))==min(min(abs(LON_BNG-XY(iii,1)))));y=y(1);
    % XY_GRID(iii,:)=[x y];
    mb=MB(x,y);
    mb=find(reg==mb);
    if strcmp(Abstraction_ID{iii,1},'A')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_GW(ano-1998,1,mb)=REG_GW(ano-1998,1,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'M')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_GW(ano-1998,2,mb)=REG_GW(ano-1998,2,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'E')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_GW(ano-1998,3,mb)=REG_GW(ano-1998,3,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'I')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_GW(ano-1998,4,mb)=REG_GW(ano-1998,4,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'P')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_GW(ano-1998,5,mb)=REG_GW(ano-1998,5,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'W')==1 && strcmp(Abstraction_ID{iii,4},'GW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_GW(ano-1998,6,mb)=REG_GW(ano-1998,6,mb)+sum(Abstraction(iii,xd));
        end
    elseif  strcmp(Abstraction_ID{iii,1},'A')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_SW(ano-1998,1,mb)=REG_SW(ano-1998,1,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'M')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_SW(ano-1998,2,mb)=REG_SW(ano-1998,2,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'E')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_SW(ano-1998,3,mb)=REG_SW(ano-1998,3,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'I')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_SW(ano-1998,4,mb)=REG_SW(ano-1998,4,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'P')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_SW(ano-1998,5,mb)=REG_SW(ano-1998,5,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'W')==1 && strcmp(Abstraction_ID{iii,4},'SW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_SW(ano-1998,6,mb)=REG_SW(ano-1998,6,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'A')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_TW(ano-1998,1,mb)=REG_TW(ano-1998,1,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'M')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_TW(ano-1998,2,mb)=REG_TW(ano-1998,2,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'E')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_TW(ano-1998,3,mb)=REG_TW(ano-1998,3,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'I')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_TW(ano-1998,4,mb)=REG_TW(ano-1998,4,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'P')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_TW(ano-1998,5,mb)=REG_TW(ano-1998,5,mb)+sum(Abstraction(iii,xd));
        end
    elseif strcmp(Abstraction_ID{iii,1},'W')==1 && strcmp(Abstraction_ID{iii,4},'TW')==1
        for ano=1999:2014
            xd=find(data(:,1)==ano);
            REG_TW(ano-1998,6,mb)=REG_TW(ano-1998,6,mb)+sum(Abstraction(iii,xd));
        end
    end
end

 
%% PLOT
nplotx=1;
nploty=2;
espvbi=0.055;
espvbs=0.045;
espvf=0.035;
esphf=0.050;
esphbe=0.05;
esphbd=0.080;
rg=[0 8];
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
c=0;

for iii=1:size(REG_GW,3)
    figure('color',[1 1 1],'position',[10 10 1300 400])
    rg=REG_GW(:,:,iii);
    subplot('position',PLOT{1,1})
    b=bar([1999:2014],rg,'stacked');
    b(1).FaceColor=[0.7 0.7 0]; % Agriculture
    b(2).FaceColor=[0.8 0.5 0.1];          % Amenity
    b(3).FaceColor=[0.1 1 0.1];          % Environmental
    b(4).FaceColor=[0.8 0.0 0.8];   % Industrial
    b(5).FaceColor=[0.7 0.2 0]; % Energy
    b(6).FaceColor=[0 0 0.9];          % Public Water Suply

    grid on
    rg=REG_SW(:,:,iii);
    subplot('position',PLOT{1,2})
    b=bar([1999:2014],rg,'stacked');
    b(1).FaceColor=[0.7 0.7 0]; % Agriculture
    b(2).FaceColor=[0.8 0.5 0.1];          % Amenity
    b(3).FaceColor=[0.1 1 0.1];          % Environmental
    b(4).FaceColor=[0.8 0 0.8];   % Industrial
    b(5).FaceColor=[0.7 0.2 0]; % Energy
    b(6).FaceColor=[0 0 0.9];          % Public Water Suply

    grid on    
    rg=REG_TW(:,:,iii);
    % subplot('position',PLOT{1,3})
    % bar([1999:2014],rg,'stacked')
    % grid on   
   % legend({'Agriculture';'Amenity';'Environmental';'Industrial, Commercial and Public Service';'Production of Energy';'Water Supply'})
    f = gcf;
    exportgraphics(f,['bargraphreg_' num2str(reg(iii)) '.png'],'Resolution',300)
    
end
