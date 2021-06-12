function fig4

fsiz = [0 0 .4 .9];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 4;
nc = 3;

h(4:6) = sim_conditioned_suppression(nr,nc,4:6,0);
h(10:12) = sim_partial_reinforcement(nr,nc,10:12,0);

for i = [1 7]
    h(i) = subplot(nr,nc,i);
    set(h(i),'visible','off');
end
h([2 3 8 9]) = [];

% --------
fs = def('fs');
fn = def('fn');
fsA = def('fsA');
xsA = -.3;def('xsA');
ysA = def('ysA');
abc = def('abc');

% abc = 'aebcdfgh';

for i= 1:length(h)
%     set(h((i)),'fontsize',fs,'fontname',fn);
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

end