%% Future Projection
clc
clear all
%% Grid
load('MB.mat')
%% Datas
c=0;
for ano=2020:2080
    for mes=1:12
        c=c+1;
        data(c,:)=[ano mes 1];
    end
end
%% Dados
ABS_GW{1,1}=zeros(size(XY,1),size(data,1));
ABS_GW{2,1}=zeros(size(XY,1),size(data,1));
ABS_GW{3,1}=zeros(size(XY,1),size(data,1));
ABS_GW{4,1}=zeros(size(XY,1),size(data,1));
ABS_GW{5,1}=zeros(size(XY,1),size(data,1));
ABS_GW{6,1}=zeros(size(XY,1),size(data,1));

ABS_GW{1,2}=zeros(size(XY,1),size(data,1));
ABS_GW{2,2}=zeros(size(XY,1),size(data,1));
ABS_GW{3,2}=zeros(size(XY,1),size(data,1));
ABS_GW{4,2}=zeros(size(XY,1),size(data,1));
ABS_GW{5,2}=zeros(size(XY,1),size(data,1));
ABS_GW{6,2}=zeros(size(XY,1),size(data,1));

ABS_GW{1,3}=zeros(size(XY,1),size(data,1));
ABS_GW{2,3}=zeros(size(XY,1),size(data,1));
ABS_GW{3,3}=zeros(size(XY,1),size(data,1));
ABS_GW{4,3}=zeros(size(XY,1),size(data,1));
ABS_GW{5,3}=zeros(size(XY,1),size(data,1));
ABS_GW{6,3}=zeros(size(XY,1),size(data,1));

ABS_SW{1,1}=zeros(size(XY,1),size(data,1));
ABS_SW{2,1}=zeros(size(XY,1),size(data,1));
ABS_SW{3,1}=zeros(size(XY,1),size(data,1));
ABS_SW{4,1}=zeros(size(XY,1),size(data,1));
ABS_SW{5,1}=zeros(size(XY,1),size(data,1));
ABS_SW{6,1}=zeros(size(XY,1),size(data,1));

ABS_SW{1,2}=zeros(size(XY,1),size(data,1));
ABS_SW{2,2}=zeros(size(XY,1),size(data,1));
ABS_SW{3,2}=zeros(size(XY,1),size(data,1));
ABS_SW{4,2}=zeros(size(XY,1),size(data,1));
ABS_SW{5,2}=zeros(size(XY,1),size(data,1));
ABS_SW{6,2}=zeros(size(XY,1),size(data,1));

ABS_SW{1,3}=zeros(size(XY,1),size(data,1));
ABS_SW{2,3}=zeros(size(XY,1),size(data,1));
ABS_SW{3,3}=zeros(size(XY,1),size(data,1));
ABS_SW{4,3}=zeros(size(XY,1),size(data,1));
ABS_SW{5,3}=zeros(size(XY,1),size(data,1));
ABS_SW{6,3}=zeros(size(XY,1),size(data,1));

ABS_TW{1,1}=zeros(size(XY,1),size(data,1));
ABS_TW{2,1}=zeros(size(XY,1),size(data,1));
ABS_TW{3,1}=zeros(size(XY,1),size(data,1));
ABS_TW{4,1}=zeros(size(XY,1),size(data,1));
ABS_TW{5,1}=zeros(size(XY,1),size(data,1));
ABS_TW{6,1}=zeros(size(XY,1),size(data,1));

ABS_TW{1,2}=zeros(size(XY,1),size(data,1));
ABS_TW{2,2}=zeros(size(XY,1),size(data,1));
ABS_TW{3,2}=zeros(size(XY,1),size(data,1));
ABS_TW{4,2}=zeros(size(XY,1),size(data,1));
ABS_TW{5,2}=zeros(size(XY,1),size(data,1));
ABS_TW{6,2}=zeros(size(XY,1),size(data,1));

ABS_TW{1,3}=zeros(size(XY,1),size(data,1));
ABS_TW{2,3}=zeros(size(XY,1),size(data,1));
ABS_TW{3,3}=zeros(size(XY,1),size(data,1));
ABS_TW{4,3}=zeros(size(XY,1),size(data,1));
ABS_TW{5,3}=zeros(size(XY,1),size(data,1));
ABS_TW{6,3}=zeros(size(XY,1),size(data,1));

%% Abertura
dd=dir('C:\Users\kl23221\OneDrive - University of Bristol\Data\Water Abstractions\Data\Future_Projection\CS-N0W_future_abstractions\*.csv');
for iii=1:size(dd,1)
    nome=['C:\Users\kl23221\OneDrive - University of Bristol\Data\Water Abstractions\Data\Future_Projection\CS-N0W_future_abstractions\' dd(iii,1).name];
    nn=split( dd(iii,1).name,'_');
    ano=str2double(nn{2,1});
    xd=find(data(:,1)==ano);
    SC=nn{3,1}(1:end-4);
    if strcmp(SC,'BaU')==1 % Scenario business as usual
        sc=1;
    elseif strcmp(SC,'Sus')==1 % Scenario ??
        sc=2;
    elseif strcmp(SC,'EG')==1 % Scenario Emisson growth 
        sc=3;
    end
    %% Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 19);

    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["OSGBEastingm", "OSGBNorthingm", "PrimaryCode", "SecondaryCode", "UseCode", "SourceCode", "Janm3month", "Febm3month", "Marm3month", "Aprm3month", "Maym3month", "Junm3month", "Julm3month", "Augm3month", "Sepm3month", "Octm3month", "Novm3month", "Decm3month", "WRZ_ID"];
    opts.VariableTypes = ["double", "double", "char", "char", "double", "char", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, ["PrimaryCode", "SecondaryCode", "SourceCode"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["PrimaryCode", "SecondaryCode", "SourceCode"], "EmptyFieldRule", "auto");

    % Import the data
    abstract = readtable(nome, opts);

    %% Convert to output type
    abstract = table2cell(abstract);
    numIdx = cellfun(@(x) ~isnan(str2double(x)), abstract);
    abstract(numIdx) = cellfun(@(x) {str2double(x)}, abstract(numIdx));
  
    %% Clear temporary variables
    
    % A - Agriculture
    % M - Amenity
    % E - Amenity
    % I - Industrial
    % P - Power
    % W - Domestic

    clear numIdx opts
    for ii=1:size(abstract,1)
        xx=XY;
        xx(:,1)=abs(xx(:,1)-abstract{ii,2});
        xx(:,2)=abs(xx(:,2)-abstract{ii,1});
        xx=sum(xx,2);
        xx=find(xx==0);
        
        if strcmp(abstract{ii,3},'A')==1 && strcmp(abstract{ii,6},'GW')==1
            ABS_GW{1,sc}(xx,xd)=ABS_GW{1,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'A')==1 && strcmp(abstract{ii,6},'SW')==1
            ABS_SW{1,sc}(xx,xd)=ABS_SW{1,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'A')==1 && strcmp(abstract{ii,6},'TW')==1
            ABS_TW{1,sc}(xx,xd)=ABS_TW{1,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'M')==1 && strcmp(abstract{ii,6},'GW')==1
            ABS_GW{2,sc}(xx,xd)=ABS_GW{2,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'M')==1 && strcmp(abstract{ii,6},'SW')==1
            ABS_SW{2,sc}(xx,xd)=ABS_SW{2,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'M')==1 && strcmp(abstract{ii,6},'TW')==1
            ABS_TW{2,sc}(xx,xd)=ABS_TW{2,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'E')==1 && strcmp(abstract{ii,6},'GW')==1
            ABS_GW{3,sc}(xx,xd)=ABS_GW{3,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'E')==1 && strcmp(abstract{ii,6},'SW')==1
            ABS_SW{3,sc}(xx,xd)=ABS_SW{3,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'E')==1 && strcmp(abstract{ii,6},'TW')==1
            ABS_TW{3,sc}(xx,xd)=ABS_TW{3,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'I')==1 && strcmp(abstract{ii,6},'GW')==1
            ABS_GW{4,sc}(xx,xd)=ABS_GW{4,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'I')==1 && strcmp(abstract{ii,6},'SW')==1
            ABS_SW{4,sc}(xx,xd)=ABS_SW{4,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'I')==1 && strcmp(abstract{ii,6},'TW')==1
            ABS_TW{4,sc}(xx,xd)=ABS_TW{4,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'P')==1 && strcmp(abstract{ii,6},'GW')==1
            ABS_GW{5,sc}(xx,xd)=ABS_GW{5,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'P')==1 && strcmp(abstract{ii,6},'SW')==1
            ABS_SW{5,sc}(xx,xd)=ABS_SW{5,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'P')==1 && strcmp(abstract{ii,6},'TW')==1
            ABS_TW{5,sc}(xx,xd)=ABS_TW{5,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'W')==1 && strcmp(abstract{ii,6},'GW')==1
            ABS_GW{6,sc}(xx,xd)=ABS_GW{6,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'W')==1 && strcmp(abstract{ii,6},'SW')==1
            ABS_SW{6,sc}(xx,xd)=ABS_SW{6,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        elseif strcmp(abstract{ii,3},'W')==1 && strcmp(abstract{ii,6},'TW')==1
            ABS_TW{6,sc}(xx,xd)=ABS_TW{6,sc}(xx,xd)+cell2mat(abstract(ii,7:18));
        end
    end
end