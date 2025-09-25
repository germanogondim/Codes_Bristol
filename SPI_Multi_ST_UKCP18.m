clc
clear all
load('UKCP18_5KM_Reg')
reg=shaperead('C:\Users\kl23221\OneDrive - University of Bristol\Data\Shapes\regions_outlines_WGS84_london_merged\regions_outlines_WGS84_london_merged.shp');
eng=[5 6 7 8 9 10 13];
area_reg=[15750 19650 8650 14975 20950 24325 15500];
% reg=reg(eng);
%% Parametros OBS
load('par_par_obs')

%% Abertura e processamento
nplotx=1;
nploty=3;
CC=1:(nplotx*nploty);
CC=reshape(CC,[nplotx nploty]);
espvbi=0.06;
espvbs=0.05;
espvf=0.055;
esphf=0.02;
esphbe=0.015;
esphbd=0.045;
[PLOT]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);

nplotx=1;
nploty=1;
CC=1:(nplotx*nploty);
CC=reshape(CC,[nplotx nploty]);
espvbi=0.075;
espvbs=0.03;
espvf=0.015;
esphf=0.00;
esphbe=0.025;
esphbd=0.035;
[PLOT2]=posicao3(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty);
load('pol_5km.mat','interp_cmap')
load('UKCP18_5KM_Reg.mat')

close all
st=[6 12 24];
for k=1:length(st)
    clear PS SPI
    data_st=P_S5KM_M_Reg{16,2};
    data_st=[data_st(:,1) data_st(:,2) ones(length(data_st(:,1)),1)];
    data_st=datenum(data_st);
    data_st=data_st(st(k):end);

    data_stv=datevec(data_st);
    p_ref=find(data_stv(:,1)>=1981 & data_stv(:,1)<=2000);
    p_nf=find(data_stv(:,1)>=2020 & data_stv(:,1)<=2040);
    p_ff=find(data_stv(:,1)>=2060 & data_stv(:,1)<=2079);

    for iii=1:size(P_S5KM_M_Reg,1)
        for ii=st(k):size(P_S5KM_M_Reg{iii,1},1)
            PS(ii-st(k)+1,:)=sum(P_S5KM_M_Reg{iii,1}(ii-st(k)+1:ii,:),1);
        end
        for ii=1:size(PS,2)
            % par=gamfit(PS(p_ref,ii));
            par=gamfit(PS(:,ii));
            % par=par_par(ii,:);
            gg=gamcdf(PS(:,ii),par(1),par(2));
            SPI(:,ii)=norminv(gg,0,1);
        end
        PS_ST{iii,1}=PS;
        UKCP18_SPI{iii,1}=SPI;
    end
    SPI_UKCP18{k,1}=UKCP18_SPI;
    SPI_UKCP18{k,2}=data_stv;
    SPI_UKCP18{k,3}=p_ref;
    SPI_UKCP18{k,4}=p_nf;
    SPI_UKCP18{k,5}=p_ff;
    c=0;
    cc=0;
    int_max_all=zeros(16*7,length(UKCP18_SPI{iii,1}(:,1)));
    int_md_all=zeros(16*7,length(UKCP18_SPI{iii,1}(:,1)));
    dm_all=zeros(16*7,length(UKCP18_SPI{iii,1}(:,1)));
    clear spi dura_hist dura_nf dura_ff nev_hist nev_nf nev_ff intmd_hist intmd_nf intmd_ff intmax_hist intmax_nf intmax_ff
     close all
    for iii=1:size(UKCP18_SPI{1,1},2)
        dm=zeros(16,length(UKCP18_SPI{iii,1}(:,1)));
        int_md=zeros(16,length(UKCP18_SPI{iii,1}(:,1)));
        int_max=zeros(16,length(UKCP18_SPI{iii,1}(:,1)));
        area=zeros(16,length(UKCP18_SPI{iii,1}(:,1)));

        md=zeros(16,length(UKCP18_SPI{iii,1}(:,1)));
        for ii=1:size(UKCP18_SPI,1)
            c=c+1;
            spi(ii,:)=UKCP18_SPI{ii,1}(:,iii);
            ss=UKCP18_SPI{ii,1}(:,iii);
            v=find(ss<0);
            F=diff(find([1 diff(v'-(1:length(v')))]));
            MM=mat2cell(v',1,[F length(v')-sum(F)]);
            for i=1:size(MM,2)
                sss=ss(MM{1,i});
                if min(sss)>-1  || length(sss)<3 || isempty(sss)==1
                    continue
                end
                cc=cc+1;
                area(ii,MM{1,i})=1;
                dm(ii,MM{1,i}(1))=length(MM{1,i});
                int_md(ii,MM{1,i}(1))=mean(sss);
                int_max(ii,MM{1,i}(1))=min(sss);

                dm_all(c,MM{1,i}(1))=length(MM{1,i});
                int_md_all(c,MM{1,i}(1))=mean(sss);
                int_max_all(c,MM{1,i}(1))=min(sss);

                md(ii,MM{1,i})=1;
                md_all(cc,:)=[length(MM{1,i}) min(sss) mean(sss) ii iii data_stv(MM{1,i}(round(length(MM{1,i})/2)),1) MM{1,i}(1) MM{1,i}(end)];
            end
            [xx,yy]=find(md==1);
            dd=dm(ii,p_ref);
            dura_hist(ii,1)=mean(dd(dd>0));
            dd=dm(ii,p_nf);
            dura_nf(ii,1)=mean(dd(dd>0));
            dd=dm(ii,p_ff);
            dura_ff(ii,1)=mean(dd(dd>0));

            dd=dm(ii,p_ref);
            nev_hist(ii,1)=length(find(dd>0));
            dd=dm(ii,p_nf);
            nev_nf(ii,1)=length(find(dd>0));
            dd=dm(ii,p_ff);
            nev_ff(ii,1)=length(find(dd>0));

            dd=int_md(ii,p_ref);
            intmd_hist(ii,1)=mean(dd(dd<0));
            dd=int_md(ii,p_nf);
            intmd_nf(ii,1)=mean(dd(dd<0));
            dd=int_md(ii,p_ff);
            intmd_ff(ii,1)=mean(dd(dd<0));

            dd=int_max(ii,p_ref);
            intmax_hist(ii,1)=mean(dd(dd<0));
            dd=int_max(ii,p_nf);
            intmax_nf(ii,1)=mean(dd(dd<0));
            dd=int_max(ii,p_ff);
            intmax_ff(ii,1)=mean(dd(dd<0));
        end
        % 
        % INTMAX=[(intmax_nf-intmax_hist)./intmax_hist (intmax_ff-intmax_hist)./intmax_hist];
        % INTMD=[(intmd_nf-intmd_hist)./intmd_hist (intmd_ff-intmd_hist)./intmd_hist];
        % DURA=[(dura_nf-dura_hist)./dura_hist (dura_ff-dura_hist)./dura_hist];
        % NEV=[(nev_nf-nev_hist)./nev_hist (nev_ff-nev_hist)./nev_hist];
        INTMAX=[(intmax_nf) (intmax_ff)];
        INTMD=[(intmd_nf) (intmd_ff)];
        DURA=[(dura_nf) (dura_ff)];
        NEV=[(nev_nf) (nev_ff)];
        INTMAX_ALL2{k,1}{iii,1}=INTMAX;
        INTMD_ALL2{k,1}{iii,1}=INTMD;
        DURA_ALL2{k,1}{iii,1}=DURA;
        NEV_ALL2{k,1}{iii,1}=NEV;
        figure('color',[1 1 1 1],'position',[100 100 1600 400])
        subplot('position',PLOT{1,1})
        % barh(flipud(NEV),'stacked')
        % set(gca,'ytick',[1:16])
        % yticklabels({'16','15','14','13','12','11','10','9','8','7','6','5','4','3','2','1'})
        boxplot(NEV,{'NF','FF'})
        title ('Number of Events')
        grid on

        subplot('position',PLOT{1,2})
        % barh(flipud(DURA),'stacked')
        % set(gca,'ytick',[1:16])
        % yticklabels({'16','15','14','13','12','11','10','9','8','7','6','5','4','3','2','1'})
        boxplot(DURA,{'NF','FF'})
        title ('Average duration')
        grid on

        subplot('position',PLOT{1,3})
        % barh(flipud(INTMAX),'stacked')
        % set(gca,'ytick',[1:16])
        % yticklabels({'16','15','14','13','12','11','10','9','8','7','6','5','4','3','2','1'})
        boxplot(INTMAX,{'NF','FF'})
        title ('Maximum average intensity')
        grid on
        % 
        % subplot('position',PLOT{1,4})
        % boxplot(INTMD,{'NF','FF'})
        % % barh(flipud(INTMAX),'stacked')
        % % legend({'NF','FF'})
        % % set(gca,'ytick',[1:16])
        % % yticklabels({'16','15','14','13','12','11','10','9','8','7','6','5','4','3','2','1'})
        % title ('Average intensity')
        grid on
        f = gcf;
        exportgraphics(f,['Caracteristicas_MD_ST' num2str(st(k),'%4.2d') reg((iii),1).geo_region '.png'],'Resolution',300)
        %% PLOT
        x1=find(spi<-3);
        x2=find(spi>-3 & spi<=-2.5);
        x3=find(spi>-2.5 & spi<=-2);
        x4=find(spi>-2 & spi<=-1.5);
        x5=find(spi>-1.5 & spi<=-1);
        x6=find(spi>-1 & spi<=-0.5);
        x7=find(spi>-0.5 & spi<=0);
        x8=find(spi>0 & spi<=0.5);
        x9=find(spi>0.5 & spi<=1);
        x10=find(spi>1 & spi<=1.5);
        x11=find(spi>1.5 & spi<=2);
        x12=find(spi>2 & spi<=2.5);
        x13=find(spi>2.5 & spi<=3);
        x14=find(spi>3);
        ss=spi;
        ss(x1)=1;
        ss(x2)=2;
        ss(x3)=3;
        ss(x4)=4;
        ss(x5)=5;
        ss(x6)=6;
        ss(x7)=7;
        ss(x8)=8;
        ss(x9)=9;
        ss(x10)=10;
        ss(x11)=11;
        ss(x12)=12;
        ss(x13)=13;
        ss(x14)=14;
        figure('color',[1 1 1 1],'position',[100 100 1600 400])
        subplot('position',PLOT2{1,1})
        imagesc(data_st,1:16,(ss))
        % imagesc(data_st,1:16,(ss),[-3 3])
        % t=pcolor(data_st,1:16,(ss));
        % t.EdgeAlpha=0.5;
        % t.LineWidth=0.1;
        colormap(interp_cmap)
        datetick('x','keeplimits')
        set(gca,'ytick',[1:16])
        f = gcf;
        exportgraphics(f,['SPI_5KM_ST' num2str(st(k),'%4.2d') reg((iii),1).geo_region '.png'],'Resolution',300)
        AREA{iii,1}=area;
    end
    INT_MAX_ALL{k,1}=int_max_all;
    INT_MD_ALL{k,1}=int_md_all;
    DM_ALL{k,1}=dm_all;
    MD_ALL{k,1}=md_all;
    DATA_STV{k,1}=data_stv;
    AREA_ST{k,1}=AREA;
end
save('RANKING_SPI_ST','INT_MAX_ALL','INT_MD_ALL','DM_ALL','DATA_STV','MD_ALL','-v7.3')
save('ranking_spi','int_max_all','int_md_all','dm_all','data_stv','md_all')
save('SPI_REG_ALL','SPI_UKCP18','-v7.3')
save('BOX_PLOT_SPATIAL','INTMAX_ALL2','INTMD_ALL2','DURA_ALL2','NEV_ALL2')

for iii=1:3
    area=0;
    for ii=1:7
    area=area+AREA_ST{iii,1}{ii,1}*area_reg(ii);
    end
    area=area/sum(area_reg);
    figure('color',[1 1 1 1],'position',[100 100 1600 400])
    subplot('position',PLOT2{1,1})
    imagesc(datenum(DATA_STV{iii,1}),1:16,(area),[0 1])
    colormap("hot")
    colorbar
    datetick('x','keeplimits')
    set(gca,'ytick',[1:16])
    f = gcf;
    exportgraphics(f,['AREA_5KM_ST' num2str(st(iii),'%4.2d') '.png'],'Resolution',300)
end

% for iii=1:size(UKCP18_SPI{iii,1},2)
%     for ii=1:size(UKCP18_SPI,1)
%         ss=UKCP18_SPI{ii,1}(:,iii);
%         v=find(ss<0);
%         F=diff(find([1 diff(v'-(1:length(v')))]));
%         MM=mat2cell(v',1,[F length(v')-sum(F)]);
%         md=zeros(size(ss,1),1)/0;
%         c=0;
%         clear md_c
%         for i=1:size(MM,2)
%             spi=ss(MM{1,i});
%             if min(spi)>-1  || length(spi)<3
%                 continue
%             end
%             c=c+1;
%             % md_c(c,:)=[data_m(min(MM{1,i}),2) data_m(min(MM{1,i}),1) data_m(max(MM{1,i}),2) data_m(max(MM{1,i}),1) length(spi) (1/normpdf(min(spi))) (1/normpdf(mean(spi)))];
%             md(MM{1,i})=1;
%         end
%         MD(:,ii)=md;
%     end
%     MD(isnan(MD))=0;
%     for ii=1:16
%         xd=find(MD(:,ii)==1);
%         for i=1:16
%         MDC(ii,i)=sum(MD(xd,i),1)/length(ss);
%         end
%     end
%     figure('color',[1 1 1 1],'position',[100 100 1600 400])
%     subplot('position',PLOT{1,1})
%     imagesc(data_st,1:16,(MDC'))
%     datetick('x','keeplimits')
%     set(gca,'ytick',[1:16])
%     f = gcf;
%     exportgraphics(f,['MDC_5KM' reg(rr(iii),1).geo_region '.png'],'Resolution',300)
%
%     MDR=zeros(length(ss),16);
%     for ii=1:16
%     xd=find(MD(:,ii)==1);
%     mm=MD(xd,:);
%     mm(:,ii)=0;
%     mm=sum(mm,2);
%     MDR(xd,ii)=mm/16;
%     end
%
%     figure('color',[1 1 1 1],'position',[100 100 1600 400])
%     subplot('position',PLOT{1,1})
%     imagesc(data_st,1:16,(MDR'),[0 0.7])
%      datetick('x','keeplimits')
%     set(gca,'ytick',[1:16])
%     colormap([[1 1 1];flipud(hot)])
%     f = gcf;
%     exportgraphics(f,['MDR_5KM' reg(rr(iii),1).geo_region '.png'],'Resolution',300)
% end
% %% Dry Days
% numColors = 256; % Número de cores na barra
% purpleColormap = [linspace(1, 0.45, numColors)', linspace(1, 0, numColors)', linspace(1, 0.45, numColors)'];
%
% for iii=1:size(P_S5KM_Reg{iii,1},2)
%     clear p0
%     for ii=1:size(P_S5KM_Reg,1)
%         p0(ii,:)=P_S5KM_M_Reg{ii,3}(:,iii);
%     end
%     x1=find(p0<5);
%     x2=find(p0>5 & p0<=10);
%     x3=find(p0>10 & p0<=15);
%     x4=find(p0>25 & p0<=20);
%     x5=find(p0>20.5 & p0<=25);
%     x6=find(p0>25 & p0<=30);
%     % x7=find(p0>-0.5 & p0<=0);
%     % x8=find(p0>0 & p0<=0.5);
%     % x9=find(p0>0.5 & p0<=1);
%     % x10=find(p0>1 & p0<=1.5);
%     % x11=find(p0>1.5 & p0<=2);
%     % x12=find(p0>2 & p0<=2.5);
%     % x13=find(p0>2.5 & p0<=3);
%     % x14=find(p0>3);
%     ss=p0;
%
%     % ss(x1)=1;
%     % ss(x2)=2;
%     % ss(x3)=3;
%     % ss(x4)=4;
%     % ss(x5)=5;
%     % ss(x6)=6;
%     % ss(x7)=7;
%     % ss(x8)=8;
%     % ss(x9)=9;
%     % ss(x10)=10;
%     % ss(x11)=11;
%     % ss(x12)=12;
%     % ss(x13)=13;
%     % ss(x14)=14;
%     figure('color',[1 1 1 1],'position',[100 100 1600 400])
%     subplot('position',PLOT{1,1})
%     imagesc(data_st,1:16,(ss),[10 31])
%
%     colormap(purpleColormap)
%     datetick('x','keeplimits')
%     set(gca,'ytick',[1:16])
%     f = gcf;
%     exportgraphics(f,['DryDays_5KM' reg(rr(iii),1).geo_region '.png'],'Resolution',300)
% end