%% SPI - Bias removed
clc
clear all

load('Fator_Correçao_Final.mat')
%% Poligonos
load('pol_5km')
[nlat,nlon]=size(LAT);
LAT_REG=ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\HadUK\data\Month\rainfall_hadukgrid_uk_1km_mon_183601-183612.nc','latitude');
LON_REG=ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\HadUK\data\Month\rainfall_hadukgrid_uk_1km_mon_183601-183612.nc','longitude');

reg=shaperead('C:\Users\kl23221\OneDrive - University of Bristol\Data\Shapes\regions_outlines_WGS84_london_merged\regions_outlines_WGS84_london_merged.shp');
ss=[5 6 7 8 9 10 13];
reg=reg(ss);
for iii=1:size(reg,1)
    in=inpolygon(LAT(1:nlat*nlon),LON(1:nlat*nlon),reg(iii,1).Y,reg(iii,1).X); % Check points inside polyon
    in=find(in);
    IN{iii,1}=in;
end
%% Arquivos
ff=[1 4 5 6 7 8 9 10 11 12 13 15 23 25 27 29];
for iii=1:length(ff)
    files_5km{iii,1}=dir(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\*.nc']);
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
for ano=1980:2080
    for mes=1:12
        % if ano==2080 && mes==12
        %     continue
        % end
        cc=cc+1;
        % data_m(cc,:)=[ano mes];
        for dia=1:30
            c=c+1;
            data_d(c,:)=[ano mes dia];
        end
    end
end
ai=[1980:10:2070];
af=[1990:10:2080];
for iii=1:12
    clc
    disp(num2str(iii*100/length(ff)))
    p_s5km=zeros(size(data_d,1),size(IN,1));

    clear ss
    for i=1:10
        pr=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\' files_5km{iii,1}(i,1).name],'pr');
        pr=rot90(pr);
        pp=zeros(size(pr,3),size(IN,1));
        xi=find(data_d(:,1)==ai(i) & data_d(:,2)==12 & data_d(:,3)==1);
        xf=find(data_d(:,1)==af(i) & data_d(:,2)==11 & data_d(:,3)==30);
        range=xi:xf;
        for ii=1:size(pr,3)
            m=data_d(range(ii),2);
            FC=zeros(size(LAT));
            FC(in_5k)=FCS_5KM{iii,1}(m,:);
            p=pr(:,:,ii).*FC;
            for k=1:size(IN,1)
                pp(ii,k)=mean(p(IN{k,1}));
            end
        end

        p_s5km(xi:xf,:)=pp;
    end
    xi=find(data_d(:,1)==1981);xi=xi(1);
    xf=find(data_d(:,1)==2079);xf=xf(end);
    p_s5km=p_s5km(xi:xf,:);
    clear PM P0;
    for ii=1:size(p_s5km,2)
        vp=p_s5km(:,ii);% Full daily time series of each grid
        vp=reshape(vp,[30 length(vp)/30]); % reorganises the time series of each grid into a matrix (30 X # months)
        v0=vp;
        v0(v0>1)=0;
        v0(v0>0)=1;
        v0=sum(v0,1);
        vpm=sum(vp,1); % Total Monthly time series
        vp=reshape(vpm,[12 length(vp)/12]); % reorganises the time series of each grid into a matrix (12 X # Years)
        clim(:,ii)=mean(vp,2);
        PM(:,ii)=vpm;
        P0(:,ii)=v0;
    end
    % PM=PM(1:end-1,:);
    % P0=P0(1:end-1,:);
    % data_m=data_m(1:end-1,:);
    P_S5KM_M_Reg{iii,1}=PM;
    P_S5KM_M_Reg{iii,2}=data_m;
    P_S5KM_Reg{iii,1}=p_s5km;
    P_S5KM_Reg{iii,2}=data_d(xi:xf,:);
    P_S5KM_M_Reg{iii,3}=P0;
end
%% 4 membros
c=0;
clear data_m
for ano=1981:2079
    for mes=1:12
        if ano==2080 && mes==12
            continue
        end
        c=c+1;
        data_m(c,:)=[ano mes];
    end
end
for iii=13:16
    clc
    disp(num2str(iii*100/length(ff)))
    data_d=[];
    %% Datas Diarias
    for i=1:10
        dd=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\' files_5km{iii,1}(i,1).name],'yyyymmdd')';
        dd=string(dd);
        datas=zeros(length(dd),3);
        for ii=1:length(dd)
            datas(ii,:)=[str2num(dd{ii,1}(1:4)) str2num(dd{ii,1}(5:6)) str2num(dd{ii,1}(7:8))];
        end
        datas=datenum(datas);
        data_d=[data_d datas'];
    end
    data_d=datevec(data_d);
    p_s5km=zeros(size(data_d,1),size(IN,1));
    clear ss
    for i=1:10
        pr=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\' files_5km{iii,1}(i,1).name],'pr');
        pr=rot90(pr);
        pp=zeros(size(pr,3),size(IN,1));
        xi=find(data_d(:,1)==ai(i) & data_d(:,2)==12 & data_d(:,3)==1);
        xf=find(data_d(:,1)==af(i) & data_d(:,2)==11 & data_d(:,3)==30);
        range=xi:xf;
        for ii=1:size(pr,3)
            m=data_d(range(ii),2);
            FC=zeros(size(LAT));
            FC(in_5k)=FCS_5KM{iii,1}(m,:);
            p=pr(:,:,ii).*FC;
            for k=1:size(IN,1)
                pp(ii,k)=mean(p(IN{k,1}));
            end
        end
 
        p_s5km(xi:xf,:)=pp;
    end
    xi=find(data_d(:,1)==1981);xi=xi(1);
    xf=find(data_d(:,1)==2079);xf=xf(end);
    p_s5km=p_s5km(xi:xf,:);
    data_d=data_d(xi:xf,:);
    c=0;
    clear PM P0
    for ano=1981:2079
        for mes=1:12
            if ano==2080 && mes==12
                continue
            end
            c=c+1;
            data_m(c,:)=[ano mes];
            xd=find(data_d(:,1)==ano & data_d(:,2)==mes);
            PM(c,:)=sum(p_s5km(xd,:),1);
            v0=p_s5km(xd,:);
            v0(v0>1)=0;
            v0(v0>0)=1;
            v0=sum(v0,1);
            P0(c,:)=v0;
        end
    end
    P_S5KM_Reg{iii,3}=P0;
    P_S5KM_M_Reg{iii,1}=PM;
    P_S5KM_M_Reg{iii,2}=data_m;
    P_S5KM_Reg{iii,1}=p_s5km;
    P_S5KM_Reg{iii,2}=data_d;
    P_S5KM_M_Reg{iii,3}=P0;
end
save('UKCP18_5KM_Reg_Bias_Corrected','P_S5KM_Reg','P_S5KM_M_Reg','-v7.3')