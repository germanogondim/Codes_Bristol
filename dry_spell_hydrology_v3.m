clc
clear all

%% arquivos
% dd=dir('*.csv');
arq(1,1).name='CAMELS_GB_hydromet_timeseries_22001_19701001-20150930.csv'
arq(2,1).name='CAMELS_GB_hydromet_timeseries_37005_19701001-20150930.csv';
arq(3,1).name='CAMELS_GB_hydromet_timeseries_39020_19701001-20150930.csv';
arq(4,1).name='CAMELS_GB_hydromet_timeseries_55026_19701001-20150930.csv';
arq(5,1).name='CAMELS_GB_hydromet_timeseries_76014_19701001-20150930.csv';
%% Plotagem
load('orange')
load('pink')
% orange=pink;
nplotx=2;
nploty=3;
espvbi=0.075;
espvbs=0.05;
espvf=0.12;
esphf=0.035;
esphbe=0.04;
esphbd=0.08;
rr=[0 8];
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
close all
%% datas
datas=datenum([1970 10 01]):datenum([2015 09 30]);
datas=datevec(datas);
c=0;
for ano=1970:2014
    for mes=1:12
        c=c+1;
        data_m(c,:)=[ano mes 1];
    end
end

%% abertura
st=[3 6 12];
for iii=1:length(st)
    figure('color',[1 1 1],'position',[100 100 1200 600])
    for k=1:5
        clear prec vazao dry_spell prec_m vazao_m aprec avazao color_ds pac qac ds_ac
        camel=camel_read(arq(k,1).name);
        prec=camel(:,1);
        vazao=camel(:,5);
        % dry spell
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
        % Escala mensal
        for ii=1:size(data_m,1)
            xd=find(datas(:,1)==data_m(ii,1) & datas(:,2)==data_m(ii,2));
            p=prec(xd);
            p=p(p>=0);
            prec_m(ii,1)=sum(p);
            v=vazao(xd);
            v=v(v>=0);
            vazao_m(ii,1)=mean(v);
            ds_m(ii,1)=sum(dry_spell(xd));
        end
        dd=data_m(st(iii):end,:);
        for ii=st(iii):length(vazao_m)
            di=data_m(ii-st(iii)+1,:);
            df=data_m(ii,:);
            xi=find(datas(:,1)==di(1) & datas(:,2)==di(2));xi=min(xi);
            xf=find(datas(:,1)==df(1) & datas(:,2)==df(2));xf=max(xf);
            qac(ii-st(iii)+1,1)=mean(vazao_m(ii-st(iii)+1:ii));
            pac(ii-st(iii)+1,1)=sum(prec_m(ii-st(iii)+1:ii));
            % ds_ac(ii-st(iii)+1,1)=sum(ds_m(ii-st(iii)+1:ii))/(st(iii)*30);
            ds_ac(ii-st(iii)+1,1)=sum(dry_spell(xi:xf))/length(xi:xf);
        end
        ds_m=ds_m(st(iii):end);
        % ds_ac=ds_ac(st(iii):end,:);
        for mes=1:12
            xd=find(dd(:,2)==mes);
            p=pac(xd);p=p(p>=0);
            v=qac(xd);v=v(v>=0);
            med_p=mean(p);dp_p=std(p);
            med_v=mean(v);dp_v=std(v);
            for ii=1:length(xd)
                aprec(xd(ii),1)=(pac(xd(ii))-med_p)/dp_p;
                avazao(xd(ii),1)=(qac(xd(ii))-med_v)/dp_v;
            end
        end

        % xs=find(dd(:,1)>2010 & dd(:,1)<=2013);
        xs=find(dd(:,1)>20);
        % xs=find(avazao<-1);


        % xs=[find(dd(:,2)==12)' find(dd(:,2)==1)' find(dd(:,2)==2)'];
        % xs2=[find(dd(:,2)==3)' find(dd(:,2)==4)' find(dd(:,2)==5)'];
        % xs3=[find(dd(:,2)==6)' find(dd(:,2)==7)' find(dd(:,2)==8)'];
        % xs4=[find(dd(:,2)==9)' find(dd(:,2)==10)' find(dd(:,2)==11)'];

        % Reclassificacao dry_spell
        % x1=find(ds_ac<0.1);
        % x2=find(ds_ac>=0.1 & ds_ac<0.20);
        % x3=find(ds_ac>=0.20 & ds_ac<0.3);
        % x4=find(ds_ac>=0.3 & ds_ac<0.4);
        % x5=find(ds_ac>=0.4 & ds_ac<0.50);
        % x6=find(ds_ac>=0.5 & ds_ac<0.60);
        % x7=find(ds_ac>=0.6 & ds_ac<0.70);
        % x8=find(ds_ac>=0.7 & ds_ac<0.80);
        % x9=find(ds_ac>=0.8);

        x1=find(ds_ac<0.05);
        x2=find(ds_ac>=0.05 & ds_ac<0.10);
        x3=find(ds_ac>=0.10 & ds_ac<0.150);
        x4=find(ds_ac>=0.15 & ds_ac<0.20);
        x5=find(ds_ac>=0.2 & ds_ac<0.250);
        x6=find(ds_ac>=0.25 & ds_ac<0.30);
        x7=find(ds_ac>=0.3 & ds_ac<0.350);
        x8=find(ds_ac>=0.35 & ds_ac<0.40);
        x9=find(ds_ac>=0.4);

        color_ds=zeros(size(ds_ac,1),3);
        color_ds(x1,:)=ones(length(x1),1)*orange(1,:);
        color_ds(x2,:)=ones(length(x2),1)*orange(2,:);
        color_ds(x3,:)=ones(length(x3),1)*orange(3,:);
        color_ds(x4,:)=ones(length(x4),1)*orange(4,:);
        color_ds(x5,:)=ones(length(x5),1)*orange(5,:);
        color_ds(x6,:)=ones(length(x6),1)*orange(6,:);
        color_ds(x7,:)=ones(length(x7),1)*orange(7,:);
        color_ds(x8,:)=ones(length(x8),1)*orange(8,:);
        color_ds(x9,:)=ones(length(x9),1)*orange(9,:);
        nome=arq(k,1).name;
        nome=split(nome,'_');
        nome=nome{5,1};
        p=aprec;p=isnan(p);
        v=avazao;v=isnan(v);
        xn=p+v;
        xn=find(xn==0);
        rr=corr(aprec(xn),avazao(xn));
        subplot('position',cell2mat(PLOT(k)))
        % plot([-4 4],[-4 4],'color',[1 0 1],'LineWidth',1.4),hold on
        % s=scatter(aprec(xs),avazao(xs),20,color_ds(xs,:),'filled');
        % text(-3,3, ['R^2 = ' num2str(rr)],'fontsize',12)
        s.MarkerEdgeColor='k';
        s.MarkerEdgeAlpha=0.1;
        for i=2:size(ds_ac,1)-1
            plot([i-1 i],[avazao(i-1) avazao(i)],'color',color_ds(i-1,:),'linewidth',1.2),hold on
        end
        % set(gca,'ylim',[-4 4],'xlim',[-4 4])
        grid on
        title(['CAMEL_ ' nome],'color','k')
        xlabel('Precipitation Anomaly')
        ylabel('Flow Anomaly')
        % ax = gca;
        % % This sets background color to black
        % ax.XColor = [1 1 1];
        % ax.YColor = [1 1 1];
        % % darkGreen = [0, 0.6, 0];
        % % ax.XColor = darkGreen;
        % ax.GridColor = 'w';
        % ax.GridAlpha = 0.5;
        % ax.FontSize = 10;
        % ax.FontWeight = 'bold';
        % hold off
        % set(gca,'Color','k')
       
    end
     f = gcf;
     exportgraphics(f,['CAMEL_flow_prec_anomaly' num2str(st(iii)) '.png'],'Resolution',300)
end