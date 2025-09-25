%% CAMEL analise teste
clc
clear all

load('camel.mat')

prec=camel(:,1);
qq=camel(:,5);
st=180;
%% dry spells
v=find(prec<=1.5);
p=zeros(size(prec));
p(v)=1;
F=diff(find([1 diff(v'-(1:length(v')))]));
dsp=mat2cell(v',1,[F length(v')-sum(F)]);
xf=find(F>=10);
dsp=cell2mat(dsp(1,xf));
dry_spell=zeros(size(prec))/0;
dry_spell(dsp)=prec(dsp);
dry_spell=dry_spell(st:end);
dry_spell(dry_spell>=0)=1; %% Identificacao dos periodos de dryspell
%% dry spells frequency
for iii=15:length(prec)
    pp=prec(iii-14:iii);
    pp(pp>1.5)=NaN;
    pp=sum(pp);
    ds(iii-14,1)=pp;
end
ds(ds>0)=1;
%% Media movel
clear ppc qqc
for iii=st:length(qq)
    ppc(iii-st+1,1)=mean(prec(iii-st+1:iii));
    qqc(iii-st+1,1)=mean(qq(iii-st+1:iii));
end
xn=find(isnan(qqc)==0);
%% Anomalia
prec=prec(st:end);
aprec=(ppc-mean(ppc))/std(ppc);
aqq=(qqc-mean(qqc(xn)))/std(qqc(xn));

v=find(aqq<-1);
p=zeros(size(aqq));
p(v)=1;
F=diff(find([1 diff(v'-(1:length(v')))]));
MM=mat2cell(v',1,[F length(v')-sum(F)]);
MM=MM(1,find(F>10));
hydro=zeros(size(aqq));
hydro(cell2mat(MM))=2;


for iii=1:size(MM,2)
    MM2{1,iii}=min(MM{1,iii})-15:max(MM{1,iii});
    hd(iii,:)=[length(MM{1,iii}) min(aqq(MM{1,iii})) length(find(dry_spell(MM{1,iii})>=0)) length(find(prec(MM{1,iii})<1)) (sum(prec(MM{1,iii})))];
end

xds=find(dry_spell==1);
for iii=1:size(MM,2)
    xxds=xds-min(MM{1,iii});
    xxds(xxds>0)=NaN;
    xxds=abs(xxds);
    xxds=xds(find(xxds==min(xxds)));
    xx=length(xxds:min(MM{1,iii}));
    x1=min(MM{1,iii})-15:min(MM{1,iii});
    x2=min(MM{1,iii})-30:min(MM{1,iii});
    hh(iii,:)=[length(MM{1,iii}) min(aqq(MM{1,iii})) sum(aqq(MM{1,iii})) length(find(dry_spell(MM{1,iii})>=0)) length(find(dry_spell(x2)>=0)) xx];
end
dry_spell(isnan(dry_spell))=0;
hydro_to_drought=hydro+dry_spell;
% hh=[hydro dry_spell]
%% 
for iii=1:size(MM,2)
    

end

