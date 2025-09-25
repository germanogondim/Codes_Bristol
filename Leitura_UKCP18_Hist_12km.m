%% Leitura UKCP18 5km
clc
clear all
%% 
clc
clear all
load('pol_12km')
% load('pol_12km')
%% Arquivos
ff=[1 4 5 6 7 8 9 10 11 12 13 15 23 25 27 29];
for iii=1:length(ff)
files_12km{iii,1}=dir(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCP18_12km\' num2str(ff(iii)) '\*.nc']);
end
%% Dates in the UKCP(12 members) format 
c=0;
for ano=1980:2020
    for mes=1:12
        for dia=1:30
            c=c+1;
            data_d(c,:)=[ano mes dia];
        end
    end
end
%% Opening 
ai=[1980 1990 2000 2010];
af=[1990 2000 2010 2020];
for iii=1:12%length(ff)
    clc
    disp(num2str(iii*100/length(ff)))
    p_s12km=zeros(size(data_d,1),length(in_12k));
    clear ss
    datas_d=[];
    for i=1:4
        pr=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCP18_12km\' num2str(ff(iii)) '\' files_12km{iii,1}(i,1).name],'pr');
        pr=rot90(pr);
        for ii=1:size(pr,3)
            p=pr(:,:,ii);
            pp(ii,:)=p(in_12k);
        end
        xi=find(data_d(:,1)==ai(i) & data_d(:,2)==12 & data_d(:,3)==1);
        xf=find(data_d(:,1)==af(i) & data_d(:,2)==11 & data_d(:,3)==30);
        p_s12km(xi:xf,:)=pp;
    end
    xi=find(data_d(:,1)==1981);xi=xi(1);
    xf=find(data_d(:,1)==2010);xf=xf(end);
    p_s12km=p_s12km(xi:xf,:);
     for ii=1:size(p_s12km,2)
            vp=p_s12km(:,ii);% Full daily time series of each grid
            vp=reshape(vp,[30 length(vp)/30]); % reorganises the time series of each grid into a matrix (30 X # months)
            vp=sum(vp,1); % Total Monthly time series 
            vp=reshape(vp,[12 length(vp)/12]); % reorganises the time series of each grid into a matrix (12 X # Years)
            clim(:,ii)=mean(vp,2);
      end
      P_S12KM_Hist{iii,1}=p_s12km;
      P_S12KM_Hist{iii,2}=data_d;
      CLIM_S12KM{iii,1}=clim;
end
%% Novos membros
for iii=13:length(ff)
    clc
    disp(num2str(iii*100/length(ff)))
    data_d=[];
    for i=1:4
        dd=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCP18_12km\' num2str(ff(iii)) '\' files_12km{iii,1}(i,1).name],'yyyymmdd')';
        dd=string(dd);
        datas=zeros(length(dd),3);
        for ii=1:length(dd)
            datas(ii,:)=[str2num(dd{ii,1}(1:4)) str2num(dd{ii,1}(5:6)) str2num(dd{ii,1}(7:8))];
        end
        datas=datenum(datas);
        data_d=[data_d datas'];
    end
    data_d=datevec(data_d);
    p_s12km=zeros(size(data_d,1),length(in_12k));
    clear ss
   
    for i=1:4
        pr=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCP18_12km\' num2str(ff(iii)) '\' files_12km{iii,1}(i,1).name],'pr');
        pr=rot90(pr);
        dd=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\UKCP18_12km\' num2str(ff(iii)) '\' files_12km{iii,1}(i,1).name],'yyyymmdd')';
        dd=string(dd);
        pp=zeros(size(pr,3),length(in_12k));
        for ii=1:size(pr,3)
            p=pr(:,:,ii);
            pp(ii,:)=p(in_12k);
        end
        xi=find(data_d(:,1)==ai(i) & data_d(:,2)==12 & data_d(:,3)==1);
        xf=find(data_d(:,1)==af(i) & data_d(:,2)==11 & data_d(:,3)==30);
        p_s12km(xi:xf,:)=pp;
    end
    xi=find(data_d(:,1)==1981);xi=xi(1);
    xf=find(data_d(:,1)==2010);xf=xf(end);
    p_s12km=p_s12km(xi:xf,:);
    data_d=data_d(xi:xf,:);
    c=0;
    for ano=1981:2010
        for mes=1:12
            c=c+1;
            xd=find(data_d(:,1)==ano & data_d(:,2)==mes);
            p_s12km_m(c,:)=sum(p_s12km(xd,:),1);
            data_m(c,:)=[ano mes];
        end
    end
    clear clim
    for mes=1:12
        xd=find(data_m(:,2)==mes);
       clim(mes,:)=mean(p_s12km_m(xd,:),1);
    end
      P_S12KM_Hist{iii,1}=p_s12km;
      P_S12KM_Hist{iii,2}=data_d;
      CLIM_S12KM{iii,1}=clim;
end
save('UKCP18_12KM','P_S12KM_Hist','CLIM_S12KM','data_d','-v7.3')

