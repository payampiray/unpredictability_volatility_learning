function fig5

fsiz = [0 0 .35 1];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 4;
nc = 2;

h(1:4) = sim_anxiety(nr,nc,1:4);
h(5:6) = sim_anxiety_piray2019(nr,nc,5:6);
h(7:8) = sim_anxiety_huang(nr,nc,7:8);

% --------
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