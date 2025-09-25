function [PLOT]=posicao(espvbi,espvbs,espvf,esphbe,esphbd,esphf,nplotx,nploty)
%%Function to help create figures with multiple plots
%nplotx - Number of rows 
%nploty - Number of coluns 
%espvbi - Bottom margin 
%espvbs - Upper margin 
%esphbe - Left margin 
%esphbd - Right margin 
%espvf  - Vertical spacing between plots
%esphf  - Horizontal spacing between plots

p1_1=esphbe;
p4=(1-espvbi -espvbs - ((nplotx-1)*espvf))/nplotx;
p3=(1-esphbe-esphbd-(((nploty*0.5)-1)*esphf))/nploty;
p2=1-espvbs-p4;
for j=1:nploty
    for i=1:nplotx
            PLOT{i,j}=[p1_1 p2 p3 p4];
            p2=p2-espvf-p4;
    end
    
    
    p1_1=p1_1+p3+esphf;
%     if j==1 || j==2 
%        p1_1=p1_1+p3+esphc;
%     elseif j==3
%        p1_1=p1_1+esphf+p3; 
%     elseif j==4|| j==5 || j==6
%       p1_1=p1_1+p3+esphc;
%     end
    p2=1-espvbs-p4;
end
end

