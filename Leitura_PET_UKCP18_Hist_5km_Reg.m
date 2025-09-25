
%% Leitura UKCP18 5km
clc
clear all
load('PET_UKCP18')
%% Poligonos
load('pol_5km')
[nlat,nlon]=size(LAT);
LAT_REG=ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\HadUK\data\Month\rainfall_hadukgrid_uk_1km_mon_183601-183612.nc','latitude');
LON_REG=ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\HadUK\data\Month\rainfall_hadukgrid_uk_1km_mon_183601-183612.nc','longitude');

reg=shaperead('C:\Users\kl23221\OneDrive - University of Bristol\Data\Shapes\regions_outlines_WGS84_london_merged\regions_outlines_WGS84_london_merged.shp');
% ss=[5 6 7 8 9 10 13];
% reg=reg(ss);
for iii=1:size(reg,1)
    in=inpolygon(LAT(1:nlat*nlon),LON(1:nlat*nlon),reg(iii,1).Y,reg(iii,1).X); % Check points inside polyon
    in=find(in);
    IN{iii,1}=in;
end
%% Data mensal
c=0;
for ano=1981:2079
    for mes=1:12
        c=c+1;
        data_m(c,:)=[ano mes];
    end
end
%% 12 Membros
c=0;
cc=0;
clear data_d
for ano=1981:2079
    for mes=1:12
        % if ano==2080 && mes==12
        %     continue
        % end
        cc=cc+1;
        data_m(cc,:)=[ano mes];
    end
end

for iii=1:12
    for ii=1:size(data_m,1)
        pet=PET_UKCP18{iii,1}(ii,:);
        MB=zeros(size(LAT_REG))/0;
        MB(in_5k)=pet;
        for i=1:size(IN,1)
        in=IN{i,1};
        pp=(MB(in));pp=mean(pp(pp>=0));
        PET(ii,i)=pp;
        end
    end
    PET_S5KM_M_Reg{iii,1}=PET;
    PET_S5KM_M_Reg{iii,2}=data_m;
end


% %% 4 membros
% c=0;
% for ano=1981:2079
%     for mes=1:12
%         if ano==2080 && mes==12
%             continue
%         end
%         c=c+1;
%         data_m(c,:)=[ano mes];
%     end
% end
% for iii=13:16
%     clc
%     disp(num2str(iii*100/length(ff)))
%     data_d=[];
%     %% Datas Diarias
%     for i=1:10
%         dd=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\' files_5km{iii,1}(i,1).name],'yyyymmdd')';
%         dd=string(dd);
%         datas=zeros(length(dd),3);
%         for ii=1:length(dd)
%             datas(ii,:)=[str2num(dd{ii,1}(1:4)) str2num(dd{ii,1}(5:6)) str2num(dd{ii,1}(7:8))];
%         end
%         datas=datenum(datas);
%         data_d=[data_d datas'];
%     end
%     data_d=datevec(data_d);
%     p_s5km=zeros(size(data_d,1),size(IN,1));
%     clear ss
%     for i=1:10
%         pr=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\' files_5km{iii,1}(i,1).name],'pr');
%         pr=rot90(pr);
%         pp=zeros(size(pr,3),size(IN,1));
%         for ii=1:size(pr,3)
%             p=pr(:,:,ii);
%             for k=1:size(IN,1)
%                 pp(ii,k)=mean(p(IN{k,1}));
%             end
%         end
%         xi=find(data_d(:,1)==ai(i) & data_d(:,2)==12 & data_d(:,3)==1);
%         xf=find(data_d(:,1)==af(i) & data_d(:,2)==11 & data_d(:,3)==30);
%         p_s5km(xi:xf,:)=pp;
%     end
%     xi=find(data_d(:,1)==1981);xi=xi(1);
%     xf=find(data_d(:,1)==2079);xf=xf(end);
%     p_s5km=p_s5km(xi:xf,:);
%     data_d=data_d(xi:xf,:);
%     c=0;
%     clear PM P0
%     for ano=1981:2079
%         for mes=1:12
%             if ano==2080 && mes==12
%                 continue
%             end
%             c=c+1;
%             data_m(c,:)=[ano mes];
%             xd=find(data_d(:,1)==ano & data_d(:,2)==mes);
%             PM(c,:)=sum(p_s5km(xd,:),1);
%             v0=p_s5km(xd,:);
%             v0(v0>1)=0;
%             v0(v0>0)=1;
%             v0=sum(v0,1);
%             P0(c,:)=v0;
%         end
%     end
%     P_S5KM_Reg{iii,3}=P0;
%     P_S5KM_M_Reg{iii,1}=PM;
%     P_S5KM_M_Reg{iii,2}=data_m;
%     P_S5KM_Reg{iii,1}=p_s5km;
%     P_S5KM_Reg{iii,2}=data_d;
%     P_S5KM_M_Reg{iii,3}=P0;
% end
save('PET_5KM_Reg','PET_S5KM_M_Reg','-v7.3')