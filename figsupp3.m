function figsupp3

fsiz = [0 0 .6 .7];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 5;
nc = 2;


h = sim_2x2_norm(nr,nc,5:10,1);
sim_2x2_norm(nr,nc,5:10,2);

% --------
fs = def('fs');
fn = def('fn');
fsA = def('fsA');
xsA = -.15;def('xsA');
ysA = def('ysA');
abc = def('abc');

for i= 1:length(h)
    set(h((i)),'fontsize',fs,'fontname',fn);
    text(xsA,ysA,abc(3+i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

end