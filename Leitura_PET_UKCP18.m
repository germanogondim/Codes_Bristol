%% Leitura PET
clc
clear all

load('pol_5km')
%% datas
c=0;
cc=0;
for ano=1981:2079
    for mes=1:12
        c=c+1;
        data_m(c,:)=[ano mes 1];
        for i=1:30
        cc=cc+1;
        data_d(cc,:)=[ano mes i];
        end
    end
end
%% Abertura
ff=[1 4 5 6 7 8 9 10 11 12 13 15];
for iii=1:12
    PP=zeros(size(data_m,1),length(in_5k));
    PP_D=zeros(size(data_d,1),length(in_5k));
    c=0;
    for ano=1981:2079
       for mes=1:12
       datai=[num2str(ano) num2str(mes,'%4.2d') '01'];
       dataf=[num2str(ano) num2str(mes,'%4.2d') '30'];

       pet=ncread(['C:\Users\kl23221\OneDrive - University of Bristol\Data\JASMIN\PET\Daily\Ens_' num2str(ff(iii),'%4.2d') '\pet_rcp85_land-cpm_uk_5km_' num2str(ff(iii),'%4.2d') '_day_' datai '-' dataf '.nc'],'pet');
       pet=rot90((pet));
       for ii=1:size(pet,3)
           p=pet(:,:,ii);
           p=p(in_5k);
           pp(ii,:)=p;
       end
       pet=sum(pp,1);
       c=c+1;
       xi=find(data_d(:,1)==ano & data_d(:,2)==mes & data_d(:,3)==1);
       xf=find(data_d(:,1)==ano & data_d(:,2)==mes & data_d(:,3)==30);
       PP(c,:)=pet;
       PET_D(xi:xf,:)=pp;
       end
    end
    PET_M_UKCP18{iii,1}=PP;
    PET_UKCP18{iii,1}=PET_D;
end
save('PET_UKCP18','PET_M_UKCP18','PET_UKCP18','-v7.3')