function [hm, hf]=plot_ci(m,e,col,alf)
if nargin<3, col = [1 .2 .2]; end
if nargin<4, alf = .3; end

if min(size(e))>1
%     the second input is errorbar
    low = e(1,:);
    high = e(2,:);
else
    low = m-e;
    high = m+e;
end

N = length(m);
if isvector(m), m=m'; end
if isvector(high), high=high'; end
if isvector(low), low=low'; end

x = 1 : N;
x2 = [x, fliplr(x)];
inBetween = [low, fliplr(high)];
hf = fill(x2, inBetween, col, 'FaceAlpha', alf, 'EdgeColor', col,'EdgeAlpha', alf); hold on;
hm = plot(x,m, 'color',col, 'LineWidth', 2);
end
