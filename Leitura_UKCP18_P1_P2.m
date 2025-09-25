%% Leitura UKCP18 5km
clc
clear all
%% 
clc
clear all
load('pol_5km')
% load('pol_12km')
%% Arquivos
ff=[1 4 5 6 7 8 9 10 11 12 13 15 23 25 27 29];
for iii=1:length(ff)
files_5km{iii,1}=dir(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\*.nc']);
end
% for iii=1:length(ff)
% files_12km{iii,1}=dir(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCP18_12km\' num2str(ff(iii)) '\*.nc']);
% end

%% P1 12 members
c=0;
clear data_d
for ano=2010:2050
    for mes=1:12
        for dia=1:30
            c=c+1;
            data_d(c,:)=[ano mes dia];
        end
    end
end
ai=[2010:10:2040];
af=[2020:10:2050];
for iii=1:12%length(ff)
    clc
    disp(num2str(iii*100/length(ff)))
    p_s5km=zeros(size(data_d,1),length(in_5k));
    clear ss
    for i=4:7
        pr=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\' files_5km{iii,1}(i,1).name],'pr');
        pr=rot90(pr);
        pp=zeros(size(pr,3),length(in_5k));
        for ii=1:size(pr,3)
            p=pr(:,:,ii);
            pp(ii,:)=p(in_5k);
        end
        xi=find(data_d(:,1)==ai(i-3) & data_d(:,2)==12 & data_d(:,3)==1);
        xf=find(data_d(:,1)==af(i-3) & data_d(:,2)==11 & data_d(:,3)==30);
        p_s5km(xi:xf,:)=pp;
    end
    xi=find(data_d(:,1)==2011);xi=xi(1);
    xf=find(data_d(:,1)==2040);xf=xf(end);
    p_s5km=p_s5km(xi:xf,:);
     for ii=1:size(p_s5km,2)
            vp=p_s5km(:,ii);% Full daily time series of each grid
            vp=reshape(vp,[30 length(vp)/30]); % reorganises the time series of each grid into a matrix (30 X # months)
            vp=sum(vp,1); % Total Monthly time series 
            vp=reshape(vp,[12 length(vp)/12]); % reorganises the time series of each grid into a matrix (12 X # Years)
            clim(:,ii)=mean(vp,2);
      end
      P_S5KM_P1{iii,1}=p_s5km;
      P_S5KM_P1{iii,2}=data_d(xi:xf,:);
      CLIM_S5KM_P1{iii,1}=clim;
end
%% P1 4 members
for iii=13:length(ff)
    clc
    disp(num2str(iii*100/length(ff)))
    data_d=[];
    for i=4:7
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
    p_s5km=zeros(size(data_d,1),length(in_5k));
    clear ss

    for i=4:7
        pr=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\' files_5km{iii,1}(i,1).name],'pr');
        pr=rot90(pr);
        pp=zeros(size(pr,3),length(in_5k));
        for ii=1:size(pr,3)
            p=pr(:,:,ii);
            pp(ii,:)=p(in_5k);
        end
        xi=find(data_d(:,1)==ai(i-3) & data_d(:,2)==12 & data_d(:,3)==1);
        xf=find(data_d(:,1)==af(i-3) & data_d(:,2)==11 & data_d(:,3)==30);
        p_s5km(xi:xf,:)=pp;
    end
    xi=find(data_d(:,1)==2011);xi=xi(1);
    xf=find(data_d(:,1)==2040);xf=xf(end);
    data_d=data_d(xi:xf,:);
    p_s5km=p_s5km(xi:xf,:);
    c=0;
    clear data_m
    for ano=2011:2040
        for mes=1:12
            c=c+1;
            xd=find(data_d(:,1)==ano & data_d(:,2)==mes);
            p_s5km_m(c,:)=sum(p_s5km(xd,:),1);
            data_m(c,:)=[ano mes];
        end
    end
    clear clim
    for mes=1:12
        xd=find(data_m(:,2)==mes);
       clim(mes,:)=mean(p_s5km_m(xd,:),1);
    end
      P_S5KM_P1{iii,1}=p_s5km;
      P_S5KM_P1{iii,2}=data_d;
      CLIM_S5KM_P1{iii,1}=clim;
end
save('UKCP18_P1','P_S5KM_P1','CLIM_S5KM_P1','data_d','-v7.3')
clear P_S5KM_P1 CLIM_S5KM_P1
%% P2 12 Members
c=0;
clear data_d
for ano=2040:2080
    for mes=1:12
        for dia=1:30
            c=c+1;
            data_d(c,:)=[ano mes dia];
        end
    end
end
ai=[2040:10:2070];
af=[2050:10:2080];
for iii=1:12%length(ff)
    clc
    disp(num2str(iii*100/length(ff)))
    p_s5km=zeros(size(data_d,1),length(in_5k));
    clear ss
    for i=7:10
        pr=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\' files_5km{iii,1}(i,1).name],'pr');
        pr=rot90(pr);
        pp=zeros(size(pr,3),length(in_5k));
        for ii=1:size(pr,3)
            p=pr(:,:,ii);
            pp(ii,:)=p(in_5k);
        end
        xi=find(data_d(:,1)==ai(i-6) & data_d(:,2)==12 & data_d(:,3)==1);
        xf=find(data_d(:,1)==af(i-6) & data_d(:,2)==11 & data_d(:,3)==30);
        p_s5km(xi:xf,:)=pp;
    end
    xi=find(data_d(:,1)==2041);xi=xi(1);
    xf=find(data_d(:,1)==2080);xf=xf(end);
    p_s5km=p_s5km(xi:xf,:);
     for ii=1:size(p_s5km,2)
            vp=p_s5km(:,ii);% Full daily time series of each grid
            vp=reshape(vp,[30 length(vp)/30]); % reorganises the time series of each grid into a matrix (30 X # months)
            vp=sum(vp,1); % Total Monthly time series 
            vp=reshape(vp,[12 length(vp)/12]); % reorganises the time series of each grid into a matrix (12 X # Years)
            clim(:,ii)=mean(vp,2);
      end
      P_S5KM_P2{iii,1}=p_s5km;
      P_S5KM_P2{iii,2}=data_d(xi:xf,:);
      CLIM_S5KM_P2{iii,1}=clim;
end
%% P2 4 Members
for iii=13:length(ff)
    clc
    disp(num2str(iii*100/length(ff)))
    data_d=[];
    for i=7:10
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
    p_s5km=zeros(size(data_d,1),length(in_5k));
    clear ss
   
    for i=7:10
        pr=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCPP_1980_2080\' num2str(ff(iii)) '\' files_5km{iii,1}(i,1).name],'pr');
        pr=rot90(pr);
        pp=zeros(size(pr,3),length(in_5k));
        for ii=1:size(pr,3)
            p=pr(:,:,ii);
            pp(ii,:)=p(in_5k);
        end
        xi=find(data_d(:,1)==ai(i-6) & data_d(:,2)==12 & data_d(:,3)==1);
        xf=find(data_d(:,1)==af(i-6) & data_d(:,2)==11 & data_d(:,3)==30);
        p_s5km(xi:xf,:)=pp;
    end
    xi=find(data_d(:,1)==2041);xi=xi(1);
    xf=find(data_d(:,1)==2080);xf=xf(end);
    p_s5km=p_s5km(xi:xf,:);
    data_d=data_d(xi:xf,:);
    c=0;
    clear data_m
    for ano=2041:2080
        for mes=1:12
            c=c+1;
            xd=find(data_d(:,1)==ano & data_d(:,2)==mes);
            p_s5km_m(c,:)=sum(p_s5km(xd,:),1);
            data_m(c,:)=[ano mes];
        end
    end
    
    clear clim
    for mes=1:12
        xd=find(data_m(:,2)==mes);
       clim(mes,:)=mean(p_s5km_m(xd,:),1);
    end
      P_S5KM_P2{iii,1}=p_s5km;
      P_S5KM_P2{iii,2}=data_d;
      CLIM_S5KM_P2{iii,1}=clim;
end
save('UKCP18_P2','P_S5KM_P2','CLIM_S5KM_P2','data_d','-v7.3')
clear P_S5KM_P2 CLIM_S5KM_P2