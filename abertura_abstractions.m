%% Analise Water Extraction
% 
clc
clear all

%% Abertura
asw=rot90(ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\Water Abstractions\Data\Gridded_Actual_Abstraction_Data_For_UoB\Gridded_Actual_Abstraction_Data_For_UoB\England_Monthly_Weighted_Abstractions_1km_Grid_SW_199901_201412.nc','abstraction'));
agw=rot90(ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\Water Abstractions\Data\Gridded_Actual_Abstraction_Data_For_UoB\Gridded_Actual_Abstraction_Data_For_UoB\England_Monthly_Weighted_Abstractions_1km_Grid_GW_199901_201412.nc','abstraction'));
%% Datas
c=0;
for ano=1999:2014
    for mes=1:12
        c=c+1;
        data(c,:)=[ano mes];
    end
end
ss=[12 1 2;3 4 5;6 7 8;9 10 11];
for iii=1:4
    xd=[];
    for ii=1:3
        xx=find(data(:,2)==ss(iii,ii));
        xd=[xd xx'];
    end
    SS{iii,1}=xd;
end
%% Poligonos
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
ASW=zeros(size(asw,3),size(reg_uk,1));
AGW=zeros(size(asw,3),size(reg_uk,1));
for iii=1:size(reg_uk,1)
    in=IN{iii,1};
    for ii=1:size(agw,3)
        aa=asw(:,:,ii);
        aa=aa(in);
        aa=sum(aa(aa>=0));
        ASW(ii,iii)=aa;
        
        aa=agw(:,:,ii);
        aa=aa(in);
        aa=sum(aa(aa>=0));
        AGW(ii,iii)=aa;
    end
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
%% Medias mensais
c=0;
clear PREC_M TEMP_M
for ano=1999:2014
    for mes=1:12
        c=c+1;
        xd=find(datas_camel(:,1)==ano & datas_camel(:,2)==mes);
        PREC_M(c,:)=sum(PREC(xd,:),1);
        TEMP_M(c,:)=mean(TEMP(xd,:),1);
        VAZAO_M(c,:)=mean(VAZAO(xd,:),1);
    end
end
%% Analise de Picos
for iii=1:size(ASW,2)
   p_min_s=islocalmin(ASW(:,iii));
   p_min_s=find(p_min_s==1);
   [~,p_max_s,~,~] = findpeaks(ASW(:,iii));

   p_min_g=islocalmin(AGW(:,iii));
   p_min_g=find(p_min_g==1);
   [~,p_max_g,~,~] = findpeaks(AGW(:,iii));
   
   max_s=zeros(size(ASW(:,iii),1),1);max_s(p_max_s)=1;
   max_g=zeros(size(ASW(:,iii),1),1);max_g(p_max_g)=3;
   min_s=zeros(size(ASW(:,iii),1),1);min_s(p_min_s)=-6;
   min_g=zeros(size(ASW(:,iii),1),1);;min_g(p_min_g)=-9;
end
%% Analise Sazonalidade
AAGW=zeros(size(AGW));
AASW=zeros(size(AGW));
ATEMP=zeros(size(AGW));
for iii=1:12
    xd=find(data(:,2)==iii);
    CLIM_AGW(iii,:)=mean(AGW(xd,:),1);
    CLIM_ASW(iii,:)=mean(ASW(xd,:),1);
    med=mean(AGW(xd,:),1);
    dp=std(AGW(xd,:),1);
    AAGW(xd,:)=(AGW(xd,:)-med)./dp;

    med=mean(ASW(xd,:),1);
    dp=std(ASW(xd,:),1);
    AASW(xd,:)=(ASW(xd,:)-med)./dp;

    med=mean(TEMP_M(xd,:),1);
    dp=std(TEMP_M(xd,:),1);
    ATEMP(xd,:)=(TEMP_M(xd,:)-med)./dp;
end
for iii=1:4
    SS_AGW(iii,:)=mean(AGW(SS{iii,1},:),1);
    SS_ASW(iii,:)=mean(ASW(SS{iii,1},:),1);
end
%% Water Abstraction year
clear YASW YAGW AYAGW AYASW AYATW
for ano=1999:2014
    xd=find(data(:,1)==ano);
    YASW(ano-1998,:)=sum(ASW(xd,:),1);
    YAGW(ano-1998,:)=sum(AGW(xd,:),1);
end
YATW=YASW+YAGW;
% Metodologia para identificar eventos de seca e como a demanda variou
% nesses periodos
for iii=1:size(YAGW,2)
    med=mean(YASW(:,iii));
    dp=std(YASW(:,iii));
    AYASW(:,iii)=(YASW(:,iii)-med)/dp;

    med=mean(YAGW(:,iii));
    dp=std(YAGW(:,iii));
    AYAGW(:,iii)=(YAGW(:,iii)-med)/dp;

    med=mean(YATW(:,iii));
    dp=std(YATW(:,iii));
    AYATW(:,iii)=(YATW(:,iii)-med)/dp;
    
end
save('Extracoes','AYATW','AYAGW','YASW','YAGW','VAZAO_M','PREC_M','TEMP_M','ASW','AGW','data','MB','AASW','AAGW','ATEMP')