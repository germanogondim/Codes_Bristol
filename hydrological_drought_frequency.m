%% Hydrological drought
clc
clear all
%% Periodo Analisado
datas=datenum([1970 10 01]):datenum([2015 09 30]);
datas=datevec(datas);
st=360;
datas_st=datas(st:end,:);
ds_ts=15;
hd_ts=-1.5;

season(1,:)=[12 1 2];
season(2,:)=[3 4 5];
season(3,:)=[6 7 8];
season(4,:)=[9 10 11];
for iii=1:4
    xs=[];
    for ii=1:3
        xd=find(datas_st(:,2)==season(iii,ii));
        xs=[xs xd'];
    end
    SS_ST{iii,1}=xs;
end
%% Analise VAZAO
arq=dir('*.csv');
HD_WPV=zeros(size(arq,1),8);
HD_WPD=zeros(size(arq,1),2);
HD_WPV_SS{1,1}=zeros(size(arq,1),8);
HD_WPV_SS{2,1}=zeros(size(arq,1),8);
HD_WPV_SS{3,1}=zeros(size(arq,1),8);
HD_WPV_SS{4,1}=zeros(size(arq,1),8);
HD_WPD_SS{1,1}=zeros(size(arq,1),4);
HD_WPD_SS{2,1}=zeros(size(arq,1),4);
for k=1:size(arq,1)
    camel=camel_read(arq(k,1).name);
    prec(:,k)=camel(:,1);
    vazao(:,k)=camel(:,5);
end

for iii=1:size(arq,1)
    vv=vazao(:,iii);
    clear v_st
    for ii=st:length(vv)
        v_st(ii-st+1,1)=mean(vv(ii-st+1:ii));
    end
    med=mean(v_st(v_st>=0));
    dp=std(v_st(v_st>=0),[],1);
    v_st=(v_st-med)/dp;
    VV_ST(:,iii)=v_st;
    % Analise Periodo Completo
    v=find((v_st)<=hd_ts);
    p=zeros(size(v_st));
    p(v)=1;
    F=diff(find([1 diff(v'-(1:length(v')))]));
    dsp=mat2cell(v',1,[F length(v')-sum(F)]);
    dsp=cell2mat(dsp);
    hd=zeros(size(v_st));
    hd(dsp)=1;
    hd_f=sum(hd)/length(hd);
    HDF(iii,1)=hd_f;
    % Season
    for ii=1:4
        ss=SS_ST{ii,1};
        HDF_SS(iii,ii)=sum(hd(ss))/length(ss);
    end
end

%% Plotagem
load('topgraphic_attributes')
n = 11;
% Inicializa o colormap
red_colormap = zeros(n, 3);
% Define a transição de preto para vermelho
for i = 1:n
    red_colormap(i, 1) = 1 -(0.5*(i-1) / (n-1)); % Componente vermelho fixo em 1
    red_colormap(i, 2) = 1 - (1.0*(i-1) / (n-1)); % Componente verde de 1 a 0
    red_colormap(i, 3) = 1 - (1.25*(i-1) / (n-1)); % Componente azul de 1 a 0
end
range_hd=0:(0.15/(size(red_colormap,1)-1)):0.15;

nplotx=1;
nploty=5;
espvbi=0.025;
espvbs=0.025;
espvf=0.035;
esphf=0.050;
esphbe=0.05;
esphbd=0.15;
rr=[0 8];
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
figure('color',[1 1 1],'position',[10 10 1200 400])
c=1;
subplot('position',PLOT{1,1})
cc=0;
clear xx yy cor_hd
for i=1:size(HDF_SS,1)
    hd=HDF(i,1);
    if isnan(hd)==1
        continue
    end
    cc=cc+1;
    xc=find(abs(range_hd-hd)==min(abs(range_hd-hd)));
    cor_hd(cc,:)=red_colormap(xc,:);
    xx(cc,:)=topgraphic_attributes(i,4);
    yy(cc,:)=topgraphic_attributes(i,3);
    % plot(topgraphic_attributes(i,4),topgraphic_attributes(i,3),'o','MarkerFaceColor',redblue(xc,:),'MarkerEdgeColor',redblue(xc,:)),hold on
end
ss=scatter(xx,yy,30,cor_hd,'filled');
ss.MarkerFaceAlpha=0.95;
ss.MarkerEdgeColor = 'k';
ss.MarkerEdgeAlpha=0.25;
axis off


for iii=1:nplotx
    for ii=1:nploty-1
        c=c+1;
        if c>5
            break
        end
        subplot('position',PLOT{iii,ii+1})
        cc=0;
        clear xx yy cor_hd
        for i=1:size(HDF_SS,1)
            hd=HDF_SS(i,c-1);
            if isnan(hd)==1
                continue
            end
            cc=cc+1;
            xc=find(abs(range_hd-hd)==min(abs(range_hd-hd)));
            cor_hd(cc,:)=red_colormap(xc,:);
            xx(cc,:)=topgraphic_attributes(i,4);
            yy(cc,:)=topgraphic_attributes(i,3);
            % plot(topgraphic_attributes(i,4),topgraphic_attributes(i,3),'o','MarkerFaceColor',redblue(xc,:),'MarkerEdgeColor',redblue(xc,:)),hold on
        end
        ss=scatter(xx,yy,30,cor_hd,'filled');
        ss.MarkerFaceAlpha=0.95;
        ss.MarkerEdgeColor = 'k';
        ss.MarkerEdgeAlpha=0.25;
        axis off
    end
end
% f = gcf;
% exportgraphics(f,['HD_WPV.png'],'Resolution',300)