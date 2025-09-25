%% Dry spells impact on hydrology
clc
clear all

%% arquivos
arq(1,1).name='CAMELS_GB_hydromet_timeseries_22001_19701001-20150930.csv';
arq(2,1).name='CAMELS_GB_hydromet_timeseries_37005_19701001-20150930.csv';
arq(3,1).name='CAMELS_GB_hydromet_timeseries_39020_19701001-20150930.csv';
arq(4,1).name='CAMELS_GB_hydromet_timeseries_55026_19701001-20150930.csv';
arq(5,1).name='CAMELS_GB_hydromet_timeseries_76014_19701001-20150930.csv';
id(1,1)=22001;
id(2,1)=37005;
id(3,1)=39020;
id(4,1)=55026;
id(5,1)=76014;
%% Atributos
load('topgraphic_attributes')
%% datas
datas=datenum([1970 10 01]):datenum([2015 09 30]);
datas=datevec(datas);
c=0;
for ano=1971:2014
    for mes=1:12
        c=c+1;
        data_m(c,:)=[ano mes 1];
    end
end
%% Weather Pattern
load('redblue')
load('orange.mat')

load('C:\Users\kl23221\OneDrive - University of Bristol\Documents\MATLAB\Weather_patterns\WP_completo')
xxi=find(datas_wp==datenum([1970 10 01]));
xxf=find(datas_wp==datenum([2015 09 30]));
wp=wp(xxi:xxf); % Recorte do WP para o periodo de estudo do CAMELS
cor_wp(1,:)=[251 180 174];
cor_wp(2,:)=[179 205 227];
cor_wp(3,:)=[204 235 197];
cor_wp(4,:)=[222 203 228];
cor_wp(5,:)=[254 217 166];
cor_wp(6,:)=[255 255 204];
cor_wp(7,:)=[229 216 189];
cor_wp(8,:)=[253 218 236];
cor_wp=cor_wp/255;
%% Abertura dos arquivos das bacias selecionadas
for k=1:5 
    camel=camel_read(arq(k,1).name);
    prec(:,k)=camel(:,1);
    vazao(:,k)=camel(:,5);
end
%% Weather Patterns
ed=1:8;
%% Media Mensais dos dados de precipitacao e vazao
for iii=1:size(data_m,1)
    xd=find(datas(:,1)==data_m(iii,1) & datas(:,2)==data_m(iii,2));
    p=prec(xd,:);
    nn=isnan(p);
    nn=sum(nn,1);
    p(isnan(p))=0;
    p=sum(p,1);
    prec_m(iii,:)=p;

    v=vazao(xd,:);
    nn=isnan(v);
    nn=sum(nn,1);
    v(isnan(v))=0;
    v=sum(v,1);
    v=v./(length(xd)-nn);
    vazao_m(iii,:)=v;
end
%% Dryspells
st=[12];
for k=1:length(st)
    for iii=1:size(prec,2)
        xd=find(topgraphic_attributes(:,1)==id(iii,1));
        area=topgraphic_attributes(xd,8);
        pp=prec(:,iii);
        v=find((pp)<=1.0);
        p=zeros(size(pp));
        p(v)=1;
        F=diff(find([1 diff(v'-(1:length(v')))]));
        dsp=mat2cell(v',1,[F length(v')-sum(F)]);
        xf=find(F>=15);
        dsp=cell2mat(dsp(1,xf));
        dd=zeros(size(pp))/0;
        dd(dsp)=prec(dsp);
        dd(dd>=0)=1; %% Identificacao dos periodos de dryspell
        dd(isnan(dd))=0;
        dry_spell(:,iii)=dd;
        vv=vazao(:,iii);
        for ii=st(k)*30:size(prec,1)
            p=prec(ii-st(k)*30+1:ii,iii);
            p=p(p>=0);
            pp_st(ii-st(k)*30+1,1)=sum(p);
            v=vazao(ii-(st(k)*30)+1:ii,1);
            v=v(v>0);
            vv_st(ii-(st(k)*30)+1,1)=mean(v);
        end
        wp_st=wp(st(k)*30:end);
        for i=30:length(wp)
            h=hist(wp(i-29:i),ed);
            h=find(h==max(h));
            h=h(1);
            wp_wp(i-29,1)=h;
        end
        prec_st=pp_st(st(k)*30:end)*(area*10^3);
        vvv_st=vv_st(st(k)*30:end)*86400;
        ds_st=dry_spell(st(k)*30:end);
        dd_st=datas(st(k)*30:end,:);
        app_st=(pp_st-mean(pp_st(pp_st>=0)))/std(pp_st(pp_st>=0));
        avv_st=(vv_st-mean(vv_st(vv_st>=0)))/std(vv_st(vv_st>=0));

        for ii=st(k):size(vazao_m,1)
            p=prec_m(ii-st(k)+1:ii,iii);
            p=p(p>=0);
            prec_st(ii-st(k)+1,:)=sum(p);

            v=vazao_m(ii-st(k)+1:ii,iii);
            v=v(v>=0);
            vazao_st(ii-st(k)+1,:)=mean(v);
        end
        aprec_st=(prec_st-mean(prec_st(prec_st>=0)))/(std(prec_st(prec_st>=0),[],1));
        avazao_st=(vazao_st-mean(vazao_st(vazao_st>=0)))/(std(vazao_st(vazao_st>=0)));
        v=find((avazao_st)<=-1);
        p=zeros(size(avazao_st,1),1);
        p(v)=1;
        F=diff(find([1 diff(v'-(1:length(v')))]));
        dsp=mat2cell(v',1,[F length(v')-sum(F)]);
        data_st=data_m(st(k):end,:);
        %% Analise de eventos
        for i=1:size(dsp,2)
            % xi=min(dsp{1,i});xi=find(dd_st(:,1)==data_st(xi,1) & dd_st(:,2)==data_st(xi,2));xi=min(xi);
            % xf=max(dsp{1,i});xf=find(dd_st(:,1)==data_st(xf,1) & dd_st(:,2)==data_st(xf,2));xf=max(xf);
            % xx=xi:xf;
            % xx=xx';
            xx=find(dd_st(:,1)>=1975 & dd_st(:,1)<=1976);
            xi=min(xx);xf=max(xf);
            xx=xx*ones(1,31);
            bb=-30:0;
            bb=ones(size(xx,1),1)*bb;
            xx=xx+bb;
            
            ds_ac=ds_st(xx);
            ds_ac=sum(ds_ac,2);
            ds_ac=ds_ac/31;
           
            xx=find(dd_st(:,1)>=1975 & dd_st(:,1)<=1976);
            xx=xx*ones(1,15);
            bb=-14:0;
            bb=ones(size(xx,1),1)*bb;
            xx=xx+bb;
            wp_wp=wp_st(xx);
            H=hist(wp_wp',ed);
            H=H-max(H,[],1);
            [xh,yh]=find(H==0);
            wp_wp=xh;

            x1=find(ds_ac==0);
            x2=find(ds_ac>0 & ds_ac<=0.10);
            x3=find(ds_ac>0.10 & ds_ac<=0.2);
            x4=find(ds_ac>0.2 & ds_ac<=0.3);
            x5=find(ds_ac>0.3 & ds_ac<=0.4);
            x6=find(ds_ac>0.4 & ds_ac<=0.5);
            x7=find(ds_ac>0.5 & ds_ac<=0.6);
            x8=find(ds_ac>0.6 & ds_ac<=0.70);
            x9=find(ds_ac>0.7 & ds_ac<=0.80);
            x10=find(ds_ac>0.8 & ds_ac<=0.9);
            x11=find(ds_ac>0.9);

            color_ds=zeros(size(ds_ac,1),3);
            color_ds(x1,:)=ones(length(x1),1)*orange2(1,:);
            color_ds(x2,:)=ones(length(x2),1)*orange2(2,:);
            color_ds(x3,:)=ones(length(x3),1)*orange2(3,:);
            color_ds(x4,:)=ones(length(x4),1)*orange2(4,:);
            color_ds(x5,:)=ones(length(x5),1)*orange2(5,:);
            color_ds(x6,:)=ones(length(x6),1)*orange2(6,:);
            color_ds(x7,:)=ones(length(x7),1)*orange2(7,:);
            color_ds(x8,:)=ones(length(x8),1)*orange2(8,:);
            color_ds(x9,:)=ones(length(x9),1)*orange2(9,:);
            color_ds(x10,:)=ones(length(x10),1)*orange2(10,:);
            color_ds(x11,:)=ones(length(x11),1)*orange2(11,:);
                   
            app=app_st(xi:xf);
            color_app=zeros(size(app,1),3);
            x1=find(app<-3);
            x2=find(app>-3 & app<-2.25);
            x3=find(app>-2.25 & app<-1.5);
            x4=find(app>-1.5 & app<-0.75);
            x5=find(app>-0.75& app<0);
            x6=find(app>0 & app<=0.75);
            x7=find(app>0.75 & app<1.5);
            x8=find(app>1.5 & app<=2.25);
            x9=find(app>2.25 & app<=3);
            x10=find(app>3);
            color_app(x1,:)=ones(length(x1),1)*redblue(1,:);
            color_app(x2,:)=ones(length(x2),1)*redblue(2,:);
            color_app(x3,:)=ones(length(x3),1)*redblue(3,:);
            color_app(x4,:)=ones(length(x4),1)*redblue(4,:);
            color_app(x5,:)=ones(length(x5),1)*redblue(5,:);
            color_app(x6,:)=ones(length(x6),1)*redblue(6,:);
            color_app(x7,:)=ones(length(x7),1)*redblue(7,:);
            color_app(x8,:)=ones(length(x8),1)*redblue(8,:);
            color_app(x9,:)=ones(length(x9),1)*redblue(9,:);
            color_app(x10,:)=ones(length(x10),1)*redblue(10,:);

            xii=xi;
            xr=min(avv_st(xi:xi+size(ds_ac,1)-1))-0.2;
            for i=1:size(ds_ac,1)-1
                plot([datenum(dd_st(xii,:)) datenum(dd_st(xii+1,:))],[avv_st(xii) avv_st(xii+1)],'color',color_ds(i,:),'linewidth',1.2),hold on
                r1=rectangle('Position',[datenum(dd_st(xii,:)) xr 1 0.1]);
                r1.FaceColor=color_app(i,:);
                r1.EdgeColor=color_app(i,:);
                r1.LineWidth=0.01;

                r2=rectangle('Position',[datenum(dd_st(xii,:)) xr-0.1 1 0.1]);
                r2.FaceColor=cor_wp(wp_wp(i),:);
                r2.EdgeColor=cor_wp(wp_wp(i),:);
                r2.LineWidth=0.01;
                xii=xii+1;
                
            end
             ax = gca;
                % This sets background color to black
                ax.XColor = [1 1 1];
                ax.YColor = [1 1 1];
                % darkGreen = [0, 0.6, 0];
                % ax.XColor = darkGreen;
                ax.GridColor = 'w';
                ax.GridAlpha = 0.5;
                ax.FontSize = 10;
                ax.FontWeight = 'bold';
                hold off
                set(gca,'Color','k')
        end
    end
end


st=[3 6 12];
for iii=1:length(st)
    for ii=st(iii):size(vazao_m,1)
        p=prec_m(ii-st(iii)+1:ii,:);
        nn=isnan(p);
        nn=sum(nn,1);
        p(isnan(p))=0;
        p=sum(p,1);
        prec_st(ii-st(iii)+1,:)=p;

        v=vazao_m(ii-st(iii)+1:ii,:);
        nn=isnan(v);
        nn=sum(nn,1);
        v(isnan(v))=0;
        v=sum(v,1);
        v=v./(st(iii)-nn);
        vazao_st(ii-st(iii)+1,:)=v;
    end
    data_st=data_m(st(iii):end,:);
    vv=vazao_st;
    nn=isnan(vv);
    nn=sum(nn,1);
    vv(isnan(vv))=0;
    vv=sum(vv,1);
    med=vv./(size(vazao_st,1)-nn);

    dp=(vazao_st-med).^2;
    nn=isnan(dp);
    nn=sum(nn,1);
    dp(isnan(dp))=0;
    dp=sum(dp,1)./(size(vazao_st,1)-nn);
    dp=dp.^0.5;
    avazao_st=(vazao_st-med)./dp;
    for ii=1:size(avazao_st,2)
        v=find((avazao_st(:,ii))<=-1);
        p=zeros(size(avazao_st,1),1);
        p(v)=1;
        F=diff(find([1 diff(v'-(1:length(v')))]));
        dsp=mat2cell(v',1,[F length(v')-sum(F)]);
        % Analise de eventos
        for i=1:size(dsp,2)
            xi=min(dsp{1,i});xi=find(datas(:,1)==data_st(xi,1) & datas(:,2)==data_st(xi,2));xi=min(xi);
            xf=max(dsp{1,i});xf=find(datas(:,1)==data_st(xf,1) & datas(:,2)==data_st(xf,2));xf=max(xf);
        end
    end
end
st=1;
iii=5;
xd=find(topgraphic_attributes(:,1)==id(iii,1));
area=topgraphic_attributes(xd,8);
TT=[prec(:,iii)*(area*10^3) vazao(:,iii)*86400];
clear TT_ST
for iii=st:size(TT,1)
    TT_ST(iii-st+1,1)=sum(TT(iii-st+1:iii,1));
    TT_ST(iii-st+1,2)=(TT(iii,2));
end
close all
plot(TT_ST(:,1)./TT_ST(:,2))

dd=pp;
dd(dd>1)=NaN;
dd(dd>=0)=1;
