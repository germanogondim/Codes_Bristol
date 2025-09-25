% %% Historical Abstractions
clc
clear all
%% Grid
load('Abstractions.mat')
xy_abs=XY;
load('MB.mat')
%% Datas
c=0;
for ano=1999:2014
    for mes=1:12
        c=c+1;
        data_h(c,:)=[ano mes 1];
    end
end
%% Abstractions
%% Dados
ABS_GW_H{1,1}=zeros(size(XY,1),size(data_h,1));
ABS_GW_H{2,1}=zeros(size(XY,1),size(data_h,1));
ABS_GW_H{3,1}=zeros(size(XY,1),size(data_h,1));
ABS_GW_H{4,1}=zeros(size(XY,1),size(data_h,1));
ABS_GW_H{5,1}=zeros(size(XY,1),size(data_h,1));
ABS_GW_H{6,1}=zeros(size(XY,1),size(data_h,1));

ABS_SW_H{1,1}=zeros(size(XY,1),size(data_h,1));
ABS_SW_H{2,1}=zeros(size(XY,1),size(data_h,1));
ABS_SW_H{3,1}=zeros(size(XY,1),size(data_h,1));
ABS_SW_H{4,1}=zeros(size(XY,1),size(data_h,1));
ABS_SW_H{5,1}=zeros(size(XY,1),size(data_h,1));
ABS_SW_H{6,1}=zeros(size(XY,1),size(data_h,1));

ABS_TW_H{1,1}=zeros(size(XY,1),size(data_h,1));
ABS_TW_H{2,1}=zeros(size(XY,1),size(data_h,1));
ABS_TW_H{3,1}=zeros(size(XY,1),size(data_h,1));
ABS_TW_H{4,1}=zeros(size(XY,1),size(data_h,1));
ABS_TW_H{5,1}=zeros(size(XY,1),size(data_h,1));
ABS_TW_H{6,1}=zeros(size(XY,1),size(data_h,1));
%%
low_flow=[70 90 100 180 200 220 230 240 250 290 320 370 430 440 450 480 630 640 650 660];
for ii=1:size(Abstraction,1)
    xx=XY; % Localizacao de todos os pontos do grid
    xx(:,1)=abs(xx(:,1)-xy_abs(ii,2));
    xx(:,2)=abs(xx(:,2)-xy_abs(ii,1));
    xx=sum(xx,2);
    xx=find(xx==0);
    aa=Abstraction(ii,:);
    ll=low_flow-Abstraction_ID{ii,3};
    ll=find(ll==0);
    if isempty(ll)==0
       aa=aa*0.3/100;
    end

    % 
    % if Abstraction_ID{ii,3}==70 || Abstraction_ID{ii,3}==90 || Abstraction_ID{ii,3}==100 || Abstraction_ID{ii,3}==180 || Abstraction_ID{ii,3}==200 ...
    %       || Abstraction_ID{ii,3}==220 || Abstraction_ID{ii,3}==230 || Abstraction_ID{ii,3}==240 || Abstraction_ID{ii,3}==250 || Abstraction_ID{ii,3}==290 || Abstraction_ID{ii,3}==320 ...
    %       || Abstraction_ID{ii,3}==370 || Abstraction_ID{ii,3}==430 || Abstraction_ID{ii,3}==440 || Abstraction_ID{ii,3}==450 || Abstraction_ID{ii,3}==480  || Abstraction_ID{ii,3}==630 ...
    %       || Abstraction_ID{ii,3}==640 || Abstraction_ID{ii,3}==650 || Abstraction_ID{ii,3}==660
    %     aa=aa*0.3/100;
    % end

    if strcmp(Abstraction_ID{ii,1},'A')==1 && strcmp(Abstraction_ID{ii,4},'GW')==1
        ABS_GW_H{1,1}(xx,:)=ABS_GW_H{1,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'A')==1 && strcmp(Abstraction_ID{ii,4},'SW')==1
        ABS_SW_H{1,1}(xx,:)=ABS_SW_H{1,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'A')==1 && strcmp(Abstraction_ID{ii,4},'TW')==1
        ABS_TW_H{1,1}(xx,:)=ABS_TW_H{1,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'M')==1 && strcmp(Abstraction_ID{ii,4},'GW')==1
        ABS_GW_H{2,1}(xx,:)=ABS_GW_H{2,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'M')==1 && strcmp(Abstraction_ID{ii,4},'SW')==1
        ABS_SW_H{2,1}(xx,:)=ABS_SW_H{2,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'M')==1 && strcmp(Abstraction_ID{ii,4},'TW')==1
        ABS_TW_H{2,1}(xx,:)=ABS_TW_H{2,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'E')==1 && strcmp(Abstraction_ID{ii,4},'GW')==1
        ABS_GW_H{3,1}(xx,:)= ABS_GW_H{3,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'E')==1 && strcmp(Abstraction_ID{ii,4},'SW')==1
        ABS_SW_H{3,1}(xx,:)=ABS_SW_H{3,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'E')==1 && strcmp(Abstraction_ID{ii,4},'TW')==1
        ABS_TW_H{3,1}(xx,:)= ABS_TW_H{3,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'I')==1 && strcmp(Abstraction_ID{ii,4},'GW')==1
        ABS_GW_H{4,1}(xx,:)= ABS_GW_H{4,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'I')==1 && strcmp(Abstraction_ID{ii,4},'SW')==1
        ABS_SW_H{4,1}(xx,:)=ABS_SW_H{4,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'I')==1 && strcmp(Abstraction_ID{ii,4},'TW')==1
        ABS_TW_H{4,1}(xx,:)=ABS_TW_H{4,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'P')==1 && strcmp(Abstraction_ID{ii,4},'GW')==1
        ABS_GW_H{5,1}(xx,:)=ABS_GW_H{5,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'P')==1 && strcmp(Abstraction_ID{ii,4},'SW')==1
        ABS_SW_H{5,1}(xx,:)=ABS_SW_H{5,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'P')==1 && strcmp(Abstraction_ID{ii,4},'TW')==1
        ABS_TW_H{5,1}(xx,:)=ABS_TW_H{5,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'W')==1 && strcmp(Abstraction_ID{ii,4},'GW')==1
        ABS_GW_H{6,1}(xx,:)=ABS_GW_H{6,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'W')==1 && strcmp(Abstraction_ID{ii,4},'SW')==1
        ABS_SW_H{6,1}(xx,:)=ABS_SW_H{6,1}(xx,:)+aa;
    elseif strcmp(Abstraction_ID{ii,1},'W')==1 && strcmp(Abstraction_ID{ii,4},'TW')==1
        ABS_TW_H{6,1}(xx,:)=ABS_TW_H{6,1}(xx,:)+aa;
    end
end

%% 
ABS_H=zeros(size(ABS_GW_H{1,1}));
for iii=1:6
    ABS_H=ABS_H+ABS_SW_H{iii,1} ;%+ ABS_GW_H{iii,1};
end

abs_nc_sw=ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\Water Abstractions\Data\Gridded_Actual_Abstraction_Data_For_UoB\Gridded_Actual_Abstraction_Data_For_UoB\England_Monthly_Weighted_Abstractions_1km_Grid_SW_199901_201412.nc','abstraction');
abs_nc_gw=ncread('C:\Users\kl23221\OneDrive - University of Bristol\Data\Water Abstractions\Data\Gridded_Actual_Abstraction_Data_For_UoB\Gridded_Actual_Abstraction_Data_For_UoB\England_Monthly_Weighted_Abstractions_1km_Grid_GW_199901_201412.nc','abstraction');
abs_nc_sw(abs_nc_sw<0)=0;
abs_nc_sw=rot90(abs_nc_sw);

abs_nc_gw(abs_nc_gw<0)=0;
abs_nc_gw=rot90(abs_nc_gw);
abs_nc=abs_nc_sw;%+abs_nc_gw;
in=find(MB>0);

for iii=1:size(abs_nc,3)
    aa=abs_nc(:,:,iii);
    aa=aa(in);
    ABS_NC(:,iii)=aa;
end

nn=ABS_NC;
nn=sum(nn,2);
nn=find(nn==0);
% ABS_NC(nn,:)=NaN;
df=abs(sum(ABS_NC,2)-sum(ABS_H,2));
df(df<0.01)=0;
% save('Historical_Abstractions','ABS_SW_H','ABS_TW_H','ABS_GW_H','data_h')


