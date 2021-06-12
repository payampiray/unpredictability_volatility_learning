function fig8

fsiz = [0 0 .7 .8];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 3;
nc = 3;

h(2:9) = sim_amygdala_costa(nr,nc,2:9,0);

for i = 1
    h(i) = subplot(nr,nc,i);
    set(h(i),'visible','off');
end

% --------
fs = def('fs');
fn = def('fn');
fsA = def('fsA');
xsA = -.15;def('xsA');
ysA = def('ysA');
abc = def('abc');

for i= 1:length(h)
    set(h((i)),'fontsize',fs,'fontname',fn);
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

end