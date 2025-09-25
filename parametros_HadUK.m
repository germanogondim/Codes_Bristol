clc
clear all
load('UKCP18_5KM_Reg')
reg=shaperead('C:\Users\kl23221\OneDrive - University of Bristol\Data\Shapes\regions_outlines_WGS84_london_merged\regions_outlines_WGS84_london_merged.shp');
ss=[5 6 7 8 9 10 13];
reg=reg(ss);
%% Parametros OBS
load('C:\Users\kl23221\OneDrive - University of Bristol\Data\HadUK\PREC_UK_Month')
load('C:\Users\kl23221\OneDrive - University of Bristol\Data\HadUK\pol_uk')
in_uk=in;
%% Grid
LAT=rot90(ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\HadUK\data\Month\rainfall_hadukgrid_uk_1km_mon_184201-184212.nc','latitude'));LAT=LAT(271:1239,200:857);nlat=size(LAT,1);
LON=rot90(ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\HadUK\data\Month\rainfall_hadukgrid_uk_1km_mon_184201-184212.nc','longitude'));LON=LON(271:1239,200:857);nlon=size(LON,2);
[nlat,nlon]=size(LAT);
%% Prec Mensal
for iii=1:size(reg,1)
    in=inpolygon(LAT(1:nlat*nlon),LON(1:nlat*nlon),reg(iii,1).Y,reg(iii,1).X); % Check points inside polyon
    in=find(in);
    IN{iii,1}=in;
end
for iii=1:size(PREC_M,2)
    MB=zeros(size(LAT))/0;
    MB(in_uk)=PREC_M(:,iii);
    for ii=1:size(IN,1)
        p=MB(IN{ii,1});
        p=p(p>=0);
        PREC_M_REG(iii,ii)=mean(p);
    end
end
st=[6 12 24];
for k=1:3
    clear P_ST_REG par_par
    for iii=st(k):size(PREC_M_REG,1)
        P_ST_REG(iii-st(k)+1,:)=sum(PREC_M_REG(iii-st(k)+1:iii,:),1);
    end
    data_st_reg=data_m(st(k):end,:);
    xx=find(data_st_reg(:,1)>=1981 & data_st_reg(:,1)<=2010);
    for iii=1:size(P_ST_REG,2)
        par_par(iii,:)=gamfit(P_ST_REG(xx,iii));
    end
    PAR_PAR{k,1}=par_par;
end
save('par_par_obs_st','PAR_PAR')