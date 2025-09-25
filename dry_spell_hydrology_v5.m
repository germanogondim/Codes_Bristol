%% Dry spells impact on hydrology
clc
clear all

%% arquivos
% arq(1,1).name='CAMELS_GB_hydromet_timeseries_22001_19701001-20150930.csv';
% arq(2,1).name='CAMELS_GB_hydromet_timeseries_37005_19701001-20150930.csv';
% arq(3,1).name='CAMELS_GB_hydromet_timeseries_39020_19701001-20150930.csv';
% arq(4,1).name='CAMELS_GB_hydromet_timeseries_55026_19701001-20150930.csv';
% arq(5,1).name='CAMELS_GB_hydromet_timeseries_76014_19701001-20150930.csv';
% arq(6,1).name='CAMELS_GB_hydromet_timeseries_8006_19701001-20150930.csv';
% arq(7,1).name='CAMELS_GB_hydromet_timeseries_50001_19701001-20150930.csv';
% arq(8,1).name='CAMELS_GB_hydromet_timeseries_28009_19701001-20150930.csv';
% arq(9,1).name='CAMELS_GB_hydromet_timeseries_39099_19701001-20150930.csv';
% arq(10,1).name='CAMELS_GB_hydromet_timeseries_39101_19701001-20150930.csv';

id(1,1)=22001;
id(2,1)=37005;
id(3,1)=39020;
id(4,1)=55026;
id(5,1)=76014;
id(6,1)=8006;
id(7,1)=50001;
id(8,1)=28009;
id(9,1)=39099;
id(10,1)=39101;

%% Atributos
load('topgraphic_attributes')
%% datas
st=360;
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
arq=dir('*.csv');

for k=1:size(arq,1)
    camel=camel_read(arq(k,1).name);
    prec(:,k)=camel(:,1);
    vazao(:,k)=camel(:,5);
end
iii=[107   109   138   193   195   219   238   289   298   310   313   329   369   445   490];
prec=prec(:,iii);
vazao=vazao(:,iii);
%% Weather Patterns
ed=1:8;
load('NAO')
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
NAO=NAO(12:end);
data_nao=data_nao(12:end,:);
xxi=find(data_nao(:,1)==1970 & data_nao(:,2)==1 & data_nao(:,3)==1); 
xxf=find(data_nao(:,1)==2015 & data_nao(:,2)==9 & data_nao(:,3)==1); 
NAO=NAO(xxi:xxf);
data_nao=data_nao(xxi:xxf,:);
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
for iii=1:size(prec,2)
    pp=prec(:,iii);
    v=find((pp)<=1.25);
    p=zeros(size(pp));
    p(v)=1;
    F=diff(find([1 diff(v'-(1:length(v')))]));
    dsp=mat2cell(v',1,[F length(v')-sum(F)]);
    dd_cum=zeros(size(pp));
    dd_vv=zeros(size(pp));
    dd_id=zeros(size(pp))/0;
    dd_ds=zeros(size(pp));
    for ii=1:size(dsp,2)
        dd_cum(dsp{1,ii})=1:length(dsp{1,ii}); 
        dd_vv(dsp{1,ii})=vazao(dsp{1,ii},iii);
        if length(dsp{1,ii})>=15
        dd_id(dsp{1,ii})=ii;
        dd_ds(dsp{1,ii})=1;
        end
    end
    dd_id=dd_id(st:end);
    dd_ds=dd_ds(st:end);
    clear vv_st
    for ii=st:length(pp)
        vv_st(ii-st+1,1)=mean(vazao(ii-st+1:ii,iii));
        pp_st(ii-st+1,1)=sum(prec(ii-st+1:ii,iii));
    end
    med=mean(vv_st(vv_st>=0));
    dp=std(vv_st(vv_st>=0));
    avv_st=(vv_st-med)/dp;
    v=find((avv_st)<=-1.0);
    p=zeros(size(pp));
    p(v)=1;
    F=diff(find([1 diff(v'-(1:length(v')))]));
    avv=mat2cell(v',1,[F length(v')-sum(F)]);
    xf=find(F>10);
    avv=avv(xf);
    pp=pp(st:end);
    clear hd
    for ii=1:size(avv,2)
        ll=(unique(dd_id(avv{1,ii})));
        ll=length(ll(ll>0));
        p=pp(avv{1,ii});p=length(find(p<=1));
        hd(ii,:)=[length(avv{1,ii}) ll sum(dd_ds(avv{1,ii})) p/length(avv{1,ii})];
    end
    v=find((avv_st)>=1.0);
    p=zeros(size(pp));
    p(v)=1;
    F=diff(find([1 diff(v'-(1:length(v')))]));
    avv=mat2cell(v',1,[F length(v')-sum(F)]);
    xf=find(F>10);
    avv=avv(xf);
    clear wet
    for ii=1:size(avv,2)
        ll=(unique(dd_id(avv{1,ii})));
        ll=length(ll(ll>0));
        p=pp(avv{1,ii});p=length(find(p<=1));
        wet(ii,:)=[length(avv{1,ii}) ll sum(dd_ds(avv{1,ii})) p/length(avv{1,ii})];
    end

    med=mean(pp_st(pp_st>=0));
    dp=std(pp_st(pp_st>=0));
    app_st=(pp_st-med)/dp;
    dd_cum=dd_cum(st:end);
    DD_CUM(:,iii)=dd_cum;
    AV_ST(:,iii)=avv_st;
    AP_ST(:,iii)=app_st;
    DD_VV(:,iii)=dd_vv;
    HD{iii,1}=hd;
    WET{iii,1}=wet;
end
%% EVENTOS
dd_st=(datas(st:end,:));
ev{1,1}=find(dd_st(:,1)>=1970 & dd_st(:,1)<=2015);
% ev{2,1}=find(dd_st(:,1)>=1990 & dd_st(:,1)<=2015);
% ev{3,1}=find(dd_st(:,1)>=1990 & dd_st(:,1)<=2000);
% ev{4,1}=find(dd_st(:,1)>=2000 & dd_st(:,1)<=2015);

dd_st=datenum(datas(st:end,:));
wp=wp(st:end,:);
for iii=20:length(wp)
    ww=wp(iii-20+1:iii);
    for ii=1:8
    wf(ii,1)=length(find(ww==ii));
    end
    wf=find(wf==max(wf));wf=wf(1);
    wp_wp(iii,1)=wf;
end   
ww=wp(1:19);
for ii=1:8
    wf(ii,1)=length(find(ww==ii));
end
wf=find(wf==max(wf));wf=wf(1);
wp_wp(1:19)=wf;
wp=wp_wp;
%% PLOTAGEM
n = 20;
% Inicializa o colormap
red_colormap = zeros(n, 3);

% Define a transição de preto para vermelho
% for i = 1:n
%     red_colormap(i, 1) = 1 -(0.5*(i-1) / (n-1)); % Componente vermelho fixo em 1
%     red_colormap(i, 2) = 1 - (1.0*(i-1) / (n-1)); % Componente verde de 1 a 0
%     red_colormap(i, 3) = 1 - (1.25*(i-1) / (n-1)); % Componente azul de 1 a 0
% end
red_colormap(1,:)=[1 1 1];
for i = 2:n
    red_colormap(i, 1) = 1 -(0.5*(i-1) / (n-1)); % Componente vermelho fixo em 1
    red_colormap(i, 2) = 0.8 - (1.0*(i-1) / (n-1)); % Componente verde de 1 a 0
    red_colormap(i, 3) = 0.8 - (1.0*(i-1) / (n-1)); % Componente azul de 1 a 0
end
red_colormap(red_colormap<0)=0;
n1=10;
clear red blue green yellow
 for i = 1:n1
    red(i, 1) = 1 -(0.7*(i-1) / (n1-1)); % Componente vermelho fixo em 1
    red(i, 2) = 1 - (1.5*(i-1) / (n1-1)); % Componente verde de 1 a 0
    red(i, 3) = 1 - (1.5*(i-1) / (n1-1)); % Componente azul de 1 a 0
 end

  for i = 1:n1
    blue(i, 3) = 1 -(0.7*(i-1) / (n1-1)); % Componente vermelho fixo em 1
    blue(i, 2) = 1 - (1.5*(i-1) / (n1-1)); % Componente verde de 1 a 0
    blue(i, 1) = 1 - (1.5*(i-1) / (n1-1)); % Componente azul de 1 a 0
  end
  for i = 1:n1
    green(i, 2) = 1 -(0.3*(i-1) / (n1-1)); % Componente vermelho fixo em 1
    green(i, 3) = 1 - (0.3*(i-1) / (n1-1)); % Componente verde de 1 a 0
    green(i, 1) = 1 - (1.5*(i-1) / (n1-1)); % Componente azul de 1 a 0
  end

    for i = 1:n1
    yellow(i, 1) = 1 -(0.3*(i-1) / (n1-1)); % Componente vermelho fixo em 1
    yellow(i, 3) = 1 - (0.3*(i-1) / (n1-1)); % Componente verde de 1 a 0
    yellow(i, 2) = 1 - (1.5*(i-1) / (n1-1)); % Componente azul de 1 a 0
  end
  greenyellow=[flipud(yellow);green];
  greenyellow(greenyellow<0)=0;
  redblue=[flipud(red);blue];
  redblue(redblue<0)=0; 
  spi_r=-2.5:(5/(size(redblue,1)-1)):2.5;
  nao_r=-3:(6/(size(greenyellow,1)-1)):3;
nplotx=1;
nploty=1;
espvbi=0.075;
espvbs=0.05;
espvf=0.05;
esphf=0.035;
esphbe=0.02;
esphbd=0.025;
rr=[0 8];
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
close all

for iii=1:size(prec,2)
    figure('color',[1 1 1],'position',[100 100 1350 400])
    for ii=1:size(ev,1)
        subplot('position',PLOT{1,ii})
        for i=1:length(ev{ii,1})-1
            rr=DD_CUM(ev{ii,1}(i),iii)+1;
            rr(rr>n)=n-1;
            if rr<5
            plot([(dd_st(ev{ii,1}(i),:)) (dd_st(ev{ii,1}(i+1),:))],[AV_ST(ev{ii,1}(i),iii) AV_ST(ev{ii,1}(i+1),iii)],'color',red_colormap(rr,:),'linewidth',0.7),hold on
            elseif rr>15
                   plot([(dd_st(ev{ii,1}(i),:)) (dd_st(ev{ii,1}(i+1),:))],[AV_ST(ev{ii,1}(i),iii) AV_ST(ev{ii,1}(i+1),iii)],'-x','MarkerSize',5,'color',red_colormap(rr,:),'linewidth',1.7),hold on
            else
                   plot([(dd_st(ev{ii,1}(i),:)) (dd_st(ev{ii,1}(i+1),:))],[AV_ST(ev{ii,1}(i),iii) AV_ST(ev{ii,1}(i+1),iii)],'color',red_colormap(rr,:),'linewidth',1.7),hold on
            end
            spi=AP_ST(ev{ii,1}(i),iii);
            xc=find(abs(spi_r-spi)==min(abs(spi_r-spi)));
            xc=xc(1);
            r1=rectangle('Position',[(dd_st(ev{ii,1}(i),:)) -3.0 1 0.5]);
            r1.FaceColor=redblue(xc,:);
            r1.EdgeColor=redblue(xc,:);
            r1.LineWidth=0.01;
            
            % r2=rectangle('Position',[(dd_st(ev{ii,1}(i),:)) -3.5 1 0.5]);
            % r2.FaceColor=cor_wp(wp(ev{ii,1}(i)),:);
            % r2.EdgeColor=cor_wp(wp(ev{ii,1}(i)),:);
            % r2.LineWidth=0.01;
            % 
            % di=datevec(dd_st(ev{ii,1}(i)));
            % xd=find(data_nao(:,1)==di(1) & data_nao(:,2)==di(2));
            % nao=NAO(xd);
            % xc=find(abs(nao_r-nao)==min(abs(nao_r-nao)));
            % xc=xc(1);
            % r3=rectangle('Position',[(dd_st(ev{ii,1}(i),:)) -3.5 1 0.5]);
            % r3.FaceColor=greenyellow(xc,:);
            % r3.EdgeColor=greenyellow(xc,:);
            % r3.LineWidth=0.01;
        end
        ax = gca;
        % This sets background color to black
        ax.XColor = [1 1 1];
        ax.YColor = [1 1 1];
        ax.XColor = [0 0 0];
        ax.YColor = [0 0 0];
        ax.GridColor = 'w';
        ax.GridAlpha = 0.5;
        ax.FontSize = 10;
        ax.FontWeight = 'bold';
        hold off
        set(gca,'Color','k','fontsize',16)
        grid on
        datetick('x','mm-yy','keeplimits')
        ylim([-4.2 4.33])
        xlim([dd_st(ev{ii,1}(1)) dd_st(ev{ii,1}(end))])
        % if ii==4
        %     colormap(red_colormap)
        %     colorbar('location','layout')
        % end
    end

    f = gcf;
    exportgraphics(f,['Figure_' num2str(id(iii)) '.png'],'Resolution',300)
end

TT=zeros(1000,30)/0;
for ii=5:30
    xd=find(DD_CUM(:,1)==ii);
    TT(1:length(xd),ii)=AV_ST(xd,1);
end


% %%
%     pp=prec(:,4);
%     v=find((pp)<=1.0);
%     p=zeros(size(pp));
%     p(v)=1;
%     F=diff(find([1 diff(v'-(1:length(v')))]));
%     xf=find(F>=5);
%     dsp=mat2cell(v',1,[F length(v')-sum(F)]);
%     for iii=5:30
%         xd=find(F==iii);
%         clear VH
%         for ii=1:length(xd)
%         VH(ii,:)=diff(vazao(dsp{1,xd(ii)},4))./vazao(dsp{1,xd(ii)}(2:end),4);
%         end
%     end
%     dsp=dsp(xf);
%     dd_cum=zeros(size(pp));
% clear dd_vv
% for ii=1:size(dsp,2)
%         dd_cum(dsp{1,ii})=1:length(dsp{1,ii});
%         dd_vv(ii,1)=max(vazao(dsp{1,ii},iii))-min(vazao(dsp{1,ii},iii));
%         dd_vv(ii,2)=(vazao(dsp{1,ii}(1),iii))-(vazao(dsp{1,ii}(end),iii));
%         dd_vv(ii,3)=length(dsp{1,ii});
%         dd_vv(ii,4)=(max(vazao(dsp{1,ii},iii))-min(vazao(dsp{1,ii},iii)))/length(dsp{1,ii});        
% end
% 
% vd=diff(vazao(:,5));
% pd=prec(2:end,5);
% v=find((vd)<=0.0);
% p=zeros(size(pp));
% p(v)=1;
% F=diff(find([1 diff(v'-(1:length(v')))]));
% dsp=mat2cell(v',1,[F length(v')-sum(F)]);
% xf=find(F>=8);
% dsp=dsp(xf);
% clear dv
% for iii=1:size(dsp,2)
%     dv(iii,:)=[abs(sum(vd(dsp{1,iii}))) sum(pd(dsp{1,iii}))];
% 
% end
% 
 
