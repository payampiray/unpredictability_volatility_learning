function figsupp6

[x1]=data_anxiety_piray2019;
x2 = huang_timeseries;

nr = 2;
nc = 1;

fs = def('fs');
fsy = def('fsy');
fn = def('fn');
fsA = def('fsA');
xsA = -.1;def('xsA');
ysA = def('ysA');
abc = def('abc');

h(1) = subplot(nr,nc,1);
plot(x1,'linewidth',2);
ylim([0 1]);
ylabel('Probability','fontsize',fsy);

h(2) = subplot(nr,nc,2);
plot(x2,'linewidth',2);
ylim([0 1]);
xlim([0 270]);
ylabel('Probability','fontsize',fsy);

for j=1:2
    text(xsA,ysA,abc(j),'fontsize',fsy,'Unit','normalized','fontname',fn,'parent',h(j));   
end

end


function [x]= huang_timeseries
N = 270;
p1 = [1 3 9 9 1 3 9 3 1 1]/13;

x = nan(N,1);
ilast = 0;
for i=1:10
    n0 = 30;
    ii = ilast + (1:n0);    
    x(ii) = p1(i);    
    ilast = ii(end);
end
x = x(1:N);

% y1 = x1 + sqrt(omega)*randn(N,1);
end
