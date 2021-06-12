function fig1

fsiz = [0 0 .4 .5];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 3;
nc = 1;

h([1 2]) = sim_feasibility(nr,nc,1:2,1);

for i = 3
    h(i) = subplot(nr,nc,i);
    set(h(i),'visible','off');
end
% h(4) = [];

% --------
fs = def('fs');
fn = def('fn');
fsA = def('fsA');
xsA = -.1;def('xsA');
ysA = def('ysA');

abc = 'abdecf';

for i= 1:length(h)
    set(h((i)),'fontsize',fs,'fontname',fn);
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

x4 = .49;
text(x4,ysA,abc(4),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(3));   

x5 = 1.1;
text(x5,ysA,abc(5),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(1));   

x6 = 1.1;
text(x6,ysA,abc(6),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(3));   


nr = 1;
nc = 2;
fsiz = [0 0 .4 .25];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

sim_feasibility(nr,nc,1:2,2);
end