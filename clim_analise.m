clc
clear all

%% Abertura
load('Historical_Abstractions.mat')
% load('Future_Abstractions.mat')
%% 
abs_gw_ht=zeros(size(ABS_GW_H{iii,1}));
abs_sw_ht=zeros(size(ABS_GW_H{iii,1}));
abs_tw_ht=zeros(size(ABS_GW_H{iii,1}));

for iii=1:6
abs_gw_ht=ABS_GW_H{iii,1}+abs_gw_ht;
abs_sw_ht=ABS_SW_H{iii,1}+abs_sw_ht;
abs_tw_ht=ABS_TW_H{iii,1}+abs_tw_ht;
end


abs_gw_hta=zeros(size(abs_gw_ht,1),length(1999:2014));
abs_sw_hta=zeros(size(abs_sw_ht,1),length(1999:2014));
abs_tw_hta=zeros(size(abs_tw_ht,1),length(1999:2014));

for ano=1999:2014
    xd=find(data_h(:,1)==ano);
    abs_gw_hta(:,ano-1998)=sum(abs_gw_ht(:,xd),2);
    abs_sw_hta(:,ano-1998)=sum(abs_sw_ht(:,xd),2);
    abs_tw_hta(:,ano-1998)=sum(abs_tw_ht(:,xd),2);
end













%% Climatologia
in=find(MB>0);
clim_abs_all=zeros(length(in),1);
for mes=1:12
    xd=find(data_h(:,2)==mes);
    xd2=find(data(:,2)==mes);
    for iii=1:6
        for i=1:size(ABS_GW{1, 1},1)
        cc=ABS_GW_H{iii,1}(i,xd);
        cc=cc(cc>0);
        clim_abs_gw_h{iii,1}(i,mes)=mean(cc);
       
        cc=ABS_SW_H{iii,1}(i,xd);
        cc=cc(cc>0);
        clim_abs_sw_h{iii,1}(i,mes)=mean(cc);
       
        cc=ABS_TW_H{iii,1}(i,xd);
        cc=cc(cc>0);
        clim_abs_tw_h{iii,1}(i,mes)=mean(cc);
        end
        c1=clim_abs_gw_h{iii,1}(:,mes);c1(isnan(c1))=0;
        c2=clim_abs_sw_h{iii,1}(:,mes);c2(isnan(c2))=0;
        c3=clim_abs_tw_h{iii,1}(:,mes);c3(isnan(c3))=0;
        clim_abs_all=clim_abs_all+c1+c2+c3;
        % for ii=1:3
        % pj_abs_gw{iii,ii}(:,xd2)=ABS_GW{iii,ii}(:,xd2)./mean(ABS_GW_H{iii,1}(:,xd),2);
        % pj_abs_sw{iii,ii}(:,xd2)=ABS_SW{iii,ii}(:,xd2)./mean(ABS_SW_H{iii,1}(:,xd),2);
        % pj_abs_tw{iii,ii}(:,xd2)=ABS_TW{iii,ii}(:,xd2)./mean(ABS_TW_H{iii,1}(:,xd),2);
        % end
    end
end
%% Limpeza
in_r=find(clim_abs_all>0);
for iii=1:6
    clim_abs_gw_h{iii,1}=clim_abs_gw_h{iii,1}(in_r,:);
    clim_abs_sw_h{iii,1}=clim_abs_sw_h{iii,1}(in_r,:);
    clim_abs_tw_h{iii,1}=clim_abs_tw_h{iii,1}(in_r,:);
    for ii=1:3
        ABS_GW{iii,ii}=ABS_GW{iii,ii}(in_r,:);
        ABS_SW{iii,ii}=ABS_SW{iii,ii}(in_r,:);
        ABS_TW{iii,ii}=ABS_TW{iii,ii}(in_r,:);
    end
end


% for iii=1:size(ABS_GW{1,1},1)
%     for mes=1:12
%         xd=find(data_h(:,2)==mes);
%         xd2=find(data(:,2)==mes);
%         pj_abs_gw{iii,}
%     end
% end
for mes=1:12
    xd=find(data_h(:,2)==mes);
    xd2=find(data(:,2)==mes);
    for iii=1:6
        for ii=1:3
        pj_abs_gw{iii,ii}(:,xd2)=(ABS_GW{iii,ii}(:,xd2)./clim_abs_gw_h{iii,1}(:,mes));
        pj_abs_sw{iii,ii}(:,xd2)=(ABS_SW{iii,ii}(:,xd2)./clim_abs_sw_h{iii,1}(:,mes));
        pj_abs_tw{iii,ii}(:,xd2)=(ABS_TW{iii,ii}(:,xd2)./ clim_abs_tw_h{iii,1}(:,mes));
        end
    end
end
