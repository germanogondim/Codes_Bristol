%% Abstractions per region
clc
clear all

%% Abertura 
load('Future_Abstractions.mat')
%% Poligonos
LAT=rot90(ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\Water Abstractions\Data\Gridded_Actual_Abstraction_Data_For_UoB\Gridded_Actual_Abstraction_Data_For_UoB\England_Monthly_Weighted_Abstractions_1km_Grid_SW_199901_201412.nc','lat'));
LON=rot90(ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\Water Abstractions\Data\Gridded_Actual_Abstraction_Data_For_UoB\Gridded_Actual_Abstraction_Data_For_UoB\England_Monthly_Weighted_Abstractions_1km_Grid_SW_199901_201412.nc','lon'));
nlat=size(LAT,1);nlon=size(LON,2);
reg_uk=shaperead('C:\Users\kl23221\OneDrive - University of Bristol\Data\Shapes\regions_outlines_WGS84_london_merged\regions_outlines_WGS84_london_merged.shp');
for iii=1:size(reg_uk,1)
    in=inpolygon(LAT(1:nlat*nlon),LON(1:nlat*nlon),reg_uk(iii,1).Y,reg_uk(iii,1).X); % Check points inside polyon
    in=find(in);
    IN{iii,1}=in;
    IN{iii,2}=reg_uk(iii,1).geo_region;
end
in=find(MB>0);
%%
for iii=1:3
    for ii=1:6
        mb_gw=zeros(size(LON));
        mb_sw=zeros(size(LON));
        mb_tw=zeros(size(LON));

        for i=1:size(data,1)
            mb_gw(in)=ABS_GW{ii,iii}(:,i);
            mb_sw(in)=ABS_SW{ii,iii}(:,i);
            mb_tw(in)=ABS_TW{ii,iii}(:,i);
            for k=1:13
            aa_gw=mb_gw(IN{k,1});
            aa_gw=sum(aa_gw(aa_gw>0));
            ABS_GW_REG{ii,iii}(k,i)=aa_gw;

            aa_sw=mb_sw(IN{k,1});
            aa_sw=sum(aa_sw(aa_sw>0));
            ABS_SW_REG{ii,iii}(k,i)=aa_sw;


            aa_tw=mb_tw(IN{k,1});
            aa_tw=sum(aa_tw(aa_sw>0));
            ABS_TW_REG{ii,iii}(k,i)=aa_tw;
            end
        end
    end
end