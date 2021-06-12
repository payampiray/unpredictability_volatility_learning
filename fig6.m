function fig6
fsiz = [.3 .3 .5 .37];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 1;
nc = 2;

h(1:2) = sim_anxiety_trait(nr,nc,1:2);

% --------
fs = def('fs');
fn = def('fn');
fsA = def('fsA');
xsA = -.15;def('xsA');
ysA = def('ysA')+0.02;
abc = def('abc');

for i= 1:length(h)
%     set(h((i)),'fontsize',fs,'fontname',fn);
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end


end