%% Abstractions Analysis
clc
clear all

%% Abertura
load('Historical_Abstractions.mat')
load('MB.mat')
load('Future_Abstractions.mat')
%% Periodo Clim
xp=find(data_h(:,1)>=2010 & data_h(:,1)<=2014);
%% Analises
pp(1,1)=78942;
pp(2,1)=154719;
pp(3,1)=92493;
pp(4,1)=95030;
pp(5,1)=76845;
pp(6,1)=152721;

%% 
gw_t=zeros(size(data_h,1),size(pp,1));
gwf_t=zeros(size(data,1),size(pp,1));
for iii=1:size(pp,1)
    gw(:,iii)=ABS_GW_H{6,1}(pp(iii,1),:);
    gwf(:,iii)=ABS_GW{6,1}(pp(iii,1),:);
    for ii=1:6
    gw_t(:,iii)=gw_t(:,iii)'+ABS_GW_H{ii,1}(pp(iii,1),:)+ABS_SW_H{ii,1}(pp(iii,1),:)+ABS_TW_H{ii,1}(pp(iii,1),:);
    gwf_t(:,iii)=gwf_t(:,iii)'+ABS_GW{ii,1}(pp(iii,1),:)+ABS_TW{ii,1}(pp(iii,1),:)+ABS_SW{ii,1}(pp(iii,1),:);
    end
end

for iii=1:12
    xc=find(data_h(:,1)>=2010 & data_h(:,2)==iii);
    for ii=1:size(pp,1)
        gg=gw(xc,ii);
        gg=gg(gg>0);
        clim_gw(iii,ii)=mean(gg);

        gg=gw_t(xc,ii);
        gg=gg(gg>0);
        clim_gwt(iii,ii)=mean(gg);
    end

end



for iii=2020:2080
    xd=find(data(:,1)==iii);
    for ii=1:size(pp,1)
    gwf_a(xd,ii)=ABS_GW{6,1}(pp(ii),xd)'./clim_gw(:,ii);
    gwft_a(xd,ii)=gwf_t(xd,(ii))./clim_gwt(:,ii);
    end
end
ABS=zeros(size(ABS_TW_H{1,1}));
for iii=1:6
    ABS=ABS+ABS_GW_H{iii,1};
end
in=find(MB>=0);
ABS_NC=zeros(size(ABS));
% nc_sw=ncread('England_Monthly_Weighted_Abstractions_1km_Grid_SW_199901_201412.nc','abstraction');nc_sw(nc_sw<0)=NaN;
nc_gw=ncread('England_Monthly_Weighted_Abstractions_1km_Grid_GW_199901_201412.nc','abstraction');nc_gw(nc_gw<0)=NaN;
abs_nc=nc_gw;
abs_nc=rot90(abs_nc);

for iii=1:size(abs_nc,3)
    aa=abs_nc(:,:,iii);
    aa=aa(in);
    ABS_NC(:,iii)=aa;
end

DF=ABS_NC-ABS;
% 
% aa(:,1)=(ABS_SW{1, 1}(136249,:));
% aa(:,2)=(ABS_SW{2, 1}(136249,:));
% aa(:,3)=(ABS_SW{3, 1}(136249,:));
% aa(:,4)=(ABS_SW{4, 1}(136249,:));
% aa(:,5)=(ABS_SW{5, 1}(136249,:));
% aa(:,6)=(ABS_SW{6, 1}(136249,:));
% aa(:,7)=(ABS_GW{1, 1}(136249,:));
% aa(:,8)=(ABS_GW{2, 1}(136249,:));
% aa(:,9)=(ABS_GW{3, 1}(136249,:));
% aa(:,10)=(ABS_GW{4, 1}(136249,:));
% aa(:,11)=(ABS_GW{5, 1}(136249,:));
% aa(:,12)=(ABS_GW{6, 1}(136249,:));
% aa(:,13)=(ABS_TW{1, 1}(136249,:));
% aa(:,14)=(ABS_TW{2, 1}(136249,:));
% aa(:,15)=(ABS_TW{3, 1}(136249,:));
% aa(:,16)=(ABS_TW{4, 1}(136249,:));
% aa(:,17)=(ABS_TW{5, 1}(136249,:));
% aa(:,18)=(ABS_TW{6, 1}(136249,:));