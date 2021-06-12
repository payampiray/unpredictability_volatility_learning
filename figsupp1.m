function figsupp1

fsiz = [0 0 .5 .5];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 2;
nc = 2;

h(1:4) = sim_volatility_task(nr,nc,1:4);


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