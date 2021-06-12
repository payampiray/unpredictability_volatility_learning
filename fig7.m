function fig7

fsiz = [0 0 .4 1];    

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 4;
nc = 2;

h(3:8) = sim_serial_prediction(nr,nc,3:8);

for i = 1
    h(i) = subplot(nr,nc,i);
    set(h(i),'visible','off');
end
h(2) = [];

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