%% HadUK 5km
clc
clear all

%%
load('pol_12km.mat')
% load('pol_12km.mat')

anoi=1981;
anof=2010;
datas_obs=datenum([anoi 01 01]):datenum([anof 12 31]);
datas_obs=datevec(datas_obs);
P_OBS12KM=zeros(size(datas_obs,1),length(in_12k));
c=0;
for ano=anoi:anof
    for mes=1:12
        c=c+1;
        xd=find(datas_obs(:,1)==ano & datas_obs(:,2)==mes);
        di=[num2str(ano) num2str(mes,'%4.2d') '01'];
        df=[num2str(ano) num2str(mes,'%4.2d') num2str(datas_obs(xd(end),3),'%4.2d')];
    
        pr=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\HadUK_12km\rainfall_hadukgrid_uk_12km_day_' di '-' df '.nc'],'rainfall');
        pr=rot90(pr);
        clear pp
        for ii=1:size(pr,3)
            p=pr(:,:,ii);
            p=p(in_12k);
            pp(ii,:)=p;
        end
        P_OBS12KM(xd,:)=pp;
        PM_OBS12KM(c,:)=sum(pp,1);
        data_m(c,:)=[ano mes];
    end
end

save('HadUK_Grid_12KM','P_OBS12KM','PM_OBS12KM','datas_obs','-v7.3')


