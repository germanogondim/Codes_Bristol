%% UKCP18 regioes
clc
clear all
%% Poligonos
load('pol_5km.mat')
[nlat,nlon]=size(LAT);
reg=shaperead('C:\Users\kl23221\OneDrive - University of Bristol\Data\Shapes\regions_outlines_WGS84_london_merged\regions_outlines_WGS84_london_merged.shp');
ss=[5 6 7 8 9 10 13];
reg=reg(ss);
for iii=1:size(reg,1)
    in=inpolygon(LAT(1:nlat*nlon),LON(1:nlat*nlon),reg(iii,1).Y,reg(iii,1).X); % Check points inside polyon
    in=find(in);
    IN{iii,1}=in;
end
%% Abertura Hist
load('UKCP18_5KM')
for iii=1:16
    data_d=P_S5KM_Hist{iii,2};
    xi=find(data_d(:,1)==1981 & data_d(:,2)==01);xi=xi(1);
    xf=find(data_d(:,1)==2010 & data_d(:,2)==12);xf=xf(end);
    c=0;
    for ano=1981:2010
        for mes=1:12
            c=c+1;
            xd=find(data_d(:,1)==ano & data_d(:,2)==mes);
            MB=zeros(size(LAT));
            MB(in_5k)=sum(P_S5KM_Hist{iii,1}(xd,:),1);
            for i=1:size(IN,1)
                p=MB(IN{i,1});
                p=p(p>=0);
                P_S5KM_REG_M{iii,1}(c,i)=mean(p);
            end
        end
    end  
end
clear P_S5KM_Hist
%% P1
load('UKCP18_P1_5KM')
for iii=1:16
    data_d=P_S5KM_P1{iii,2};
    c=size(P_S5KM_REG_M{iii,1},1);
    for ano=2012:2042
        for mes=1:12
            c=c+1;
            xd=find(data_d(:,1)==ano & data_d(:,2)==mes);
            MB=zeros(size(LAT));
            MB(in_5k)=sum(P_S5KM_P1{iii,1}(xd,:),1);
            for i=1:size(IN,1)
                p=MB(IN{i,1});
                p=p(p>=0);
                P_S5KM_REG_M{iii,1}(c,i)=mean(p);
            end
        end
    end  
end
clear P_S5KM_P1
%% P2
load('UKCP18_P1_5KM')
for iii=1:16
    data_d=P_S5KM_P1{iii,2};
    c=size(P_S5KM_REG_M{iii,1},1);
    for ano=2012:2042
        for mes=1:12
            c=c+1;
            xd=find(data_d(:,1)==ano & data_d(:,2)==mes);
            MB=zeros(size(LAT));
            MB(in_5k)=sum(P_S5KM_P2{iii,1}(xd,:),1);
            for i=1:size(IN,1)
                p=MB(IN{i,1});
                p=p(p>=0);
                P_S5KM_REG_M{iii,1}(c,i)=mean(p);
            end
        end
    end  
end
save('UKCP18_5KM_Reg','P_S5KM_REG_M')