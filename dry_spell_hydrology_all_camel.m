%% Dry_Spell Hydrology
% 1 - Identificar os extremos negativos hidrologicos
% 2 - Atribuir a ocorrencia dos extremos hidrologicos negativos (EHN) a dry spells (DS)
% 3 - Avaliar a relacao entre EHN e as caracteristicas dos DS

clc
clear all
%% arquivos
% dd=dir('*.csv');
arq(1,1).name='CAMELS_GB_hydromet_timeseries_22001_19701001-20150930.csv';
arq(2,1).name='CAMELS_GB_hydromet_timeseries_37005_19701001-20150930.csv';
arq(3,1).name='CAMELS_GB_hydromet_timeseries_39020_19701001-20150930.csv';
arq(4,1).name='CAMELS_GB_hydromet_timeseries_55026_19701001-20150930.csv';
arq(5,1).name='CAMELS_GB_hydromet_timeseries_76014_19701001-20150930.csv';

c=0;
%% Plotagemmm=zeros(size(MB))/0;
nplotx=2;
nploty=3;
espvbi=0.035;
espvbs=0.07;
espvf=0.08;
esphf=0.035;
esphbe=0.04;
esphbd=0.08;
rr=[0 8];
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
close all
figure('color',[1 1 1 1],'position',[100 100 1200 600])
%% Abertura aquivos
for k=1:size(arq,1)
    clear avazao aprec ppc vazao prec qqc
    camel=camel_read(arq(k,1).name);
    prec=camel(:,1);
    vazao=camel(:,5);
    datas=datenum([1970 10 01]):datenum([2015 09 30]);
    datas=datevec(datas);
    ss(1,:)=[12 1 2];
    ss(2,:)=[3 4 5];
    ss(3,:)=[6 7 8];
    ss(4,:)=[9 10 11];
    %% Identificacao dos ENH e dry spells
    st=360;
    datas=datas(st:end,:);
    for iii=st:length(vazao)
        ppc(iii-st+1,1)=sum(prec(iii-st+1:iii));
        qqc(iii-st+1,1)=mean(vazao(iii-st+1:iii));
    end
    xn=find(isnan(qqc)==0);
    prec=prec(st:end);
    vazao=vazao(st:end);

    for iii=1:4
        xd=[];
        for ii=1:3
            xd=[xd find(datas(:,2)==ss(iii,ii))'];
        end
        SS{iii,1}=xd;
    end
    for iii=1:4
        xs=SS{iii,1};
        vv=qqc(xs);
        vv=vv(vv>=0);
        med=mean(vv);
        dp=std(vv);
        avazao(xs)=(qqc(xs)-med)/dp;
    end
    avazao=avazao';

    % for mes=1:12
    %     xd=find(datas(:,2)==mes);
    %     vv=qqc(xd);
    %     vv=vv(vv>=0);
    %     med=mean(vv);
    %     dp=std(vv);
    %     for ii=1:length(xd)
    %         avazao(xd(ii),1)=(qqc(xd(ii))-med)/dp;
    %     end
    %
    % end
    aprec=(ppc-mean(ppc))/std(ppc);
    avazao2=(qqc-mean(qqc(xn)))/std(qqc(xn),[],1);

    v=find(avazao<-1);
    p=zeros(size(avazao));
    p(v)=1;
    F=diff(find([1 diff(v'-(1:length(v')))]));
    MM=mat2cell(v',1,[F length(v')-sum(F)]);
    MM=MM(1,find(F>=15));
    ENH=zeros(size(avazao));
    ENH(cell2mat(MM))=2;
    xx=find(ENH==2);
    qref=zeros(size(vazao))/0;
    qref(xx)=vazao(xx);


    v=find((ppc/st)<=1.5);
    p=zeros(size(ppc));
    p(v)=1;
    F=diff(find([1 diff(v'-(1:length(v')))]));
    dsp=mat2cell(v',1,[F length(v')-sum(F)]);
    xf=find(F>=15);
    dsp=cell2mat(dsp(1,xf));
    dry_spell_st=zeros(size(ppc))/0;
    dry_spell_st(dsp)=ppc(dsp)/st;
    dry_spell_st(dry_spell_st>=0)=1; %% Identificacao dos periodos de dryspell
    dry_spell_st(isnan(dry_spell_st))=0;

    v=find((prec)<=1.5);
    p=zeros(size(prec));
    p(v)=1;
    F=diff(find([1 diff(v'-(1:length(v')))]));
    dsp=mat2cell(v',1,[F length(v')-sum(F)]);
    xf=find(F>=15);
    dsp=cell2mat(dsp(1,xf));
    dry_spell=zeros(size(prec))/0;
    dry_spell(dsp)=prec(dsp);
    dry_spell(dry_spell>=0)=1; %% Identificacao dos periodos de dryspell
    dry_spell(isnan(dry_spell))=0;
    %% Ref Prec
    for iii=15:length(prec)
        p15(iii-14,1)=sum(prec(iii-14:iii));
    end
    p15=mean(p15);
    for iii=30:length(prec)
        p30(iii-29,1)=sum(prec(iii-29:iii));
    end
    p30=mean(p30);
    %% Analise dos ENH
    for iii=1:size(MM,2)
        xi=min(MM{1,iii});
        if xi<30
            continue
        end
        c=c+1;
        xf=max(MM{1,iii});
        dura=length(MM{1,iii});
        int=sum(avazao(MM{1,iii}));
        ds=sum(dry_spell_st(MM{1,iii}));
        ds15=sum(dry_spell_st(xi-15:xf));
        ds30=sum(dry_spell_st(xi-30:xf));
        
        pr=sum(prec(MM{1,iii}));
        pr15=sum(prec(xi-15:xf));
        pr30=sum(prec(xi-30:xf));
        dd=prec(xi:xf);
        dd=length(find(dd<=1));
        dd15=prec(xi-15:xf);
        dd15=length(find(dd15<=1.0));
        dd30=prec(xi-30:xf);
        dd30=length(find(dd30<=1.0));
        HD(c,:)=[dura int ds ds15 ds30 pr pr15 pr30 dd dd15 dd30];
    end
    nome=arq(k,1).name;
    nome=split(nome,'_');
    nome=nome{5,1};
    subplot('position',cell2mat(PLOT(k)))
    plot(datenum(datas),avazao)
    dateaxis()
    ylim([-2.5 4.5])
    grid on
    title(['CAMEL_ ' nome])   
end
f = gcf;  
exportgraphics(f,['CAMEL_flowanomaly' num2str(st) '.png'],'Resolution',300)
%% Plotagemmm=zeros(size(MB))/0;
nplotx=3;
nploty=3;
espvbi=0.025;
espvbs=0.09;
espvf=0.08;
esphf=0.05;
esphbe=0.085;
esphbd=0.1;
rr=[0 8];
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
xx=find(HD(:,4)>0);
figure('color',[1 1 1 1],'position',[100 100 800 1200])
subplot('position',PLOT{1,1})
scatter(HD(xx,1),HD(xx,3))
ylim([0 400])
grid on
subplot('position',PLOT{1,2})
scatter(HD(xx,1),HD(xx,4))
ylim([0 400])
grid on
subplot('position',PLOT{1,3})
scatter(HD(xx,1),HD(xx,5))
ylim([0 400])
grid on

subplot('position',PLOT{2,1})
scatter(HD(xx,1),HD(xx,6))
grid on
subplot('position',PLOT{2,2})
scatter(HD(xx,1),HD(xx,7))
grid on
subplot('position',PLOT{2,3})
scatter(HD(xx,1),HD(xx,8))
grid on

subplot('position',PLOT{3,1})
scatter(HD(xx,1),HD(xx,9))
ylim([0 300])
grid on
subplot('position',PLOT{3,2})
scatter(HD(xx,1),HD(xx,10))
grid on
subplot('position',PLOT{3,3})
scatter(HD(xx,1),HD(xx,11))
grid on


f = gcf;
exportgraphics(f,['Figure_2_st' num2str(st) '.png'],'Resolution',300)