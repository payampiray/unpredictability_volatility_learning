function hb = errorbarKxN(mx,ex,labels,colmap,barwidth,ylbl)
% examples
% 
% % factorial
% close all;
% figure;
% mx = [0.4 0.8 1.2; 0.3 0.6 .9];
% ex = [0.15 0.1 .02; 0.08 0.1 .1];
% labels = {'group 1','group 2','group 3'};
% colmap = [0.2000    0.6000    1.0000; 0.0500    0.1500    0.3000];
% barwidth = .2;
% run(mx,ex,labels,colmap,barwidth);
% 
% % factorial 2x2
% m = [0.5191    0.7112    0.3213    0.4359];
% % factors: NlVl, NlVh, NhVl, NhVh
% mx = [m([1 2])' m([3 4])'];
% ex = mx*0;
% labels = {'low obs noise','high obs noise'};
% figure;
% errorbarKxN(mx,ex,labels);
% legend({'low volatility','high volatility'},'fontsize',15)
% 
% % a 3-column bar with different colors
% figure;
% mx = [0.4 0.8 1.2]';
% ex = [0.15 0.1 .02]';
% labels = {'factor 1','factor 2','factor 3'};
% colmap = [0.2000    0.6000    1.0000; 0.0500    0.1500    0.3000; 0.5500    0.1500    0.3000];
% barwidth = .05;
% run(mx,ex,labels,colmap,barwidth);
% 
% % a 3-column bar with different colors
% figure;
% mx = [0.4 0.8 1.2];
% ex = [0.15 0.1 .02];
% colmap = [0.2000    0.6000    1.0000];
% barwidth = .2;
% labels = {'group 1','group 2','group 3'};
% run(mx,ex,labels,colmap,barwidth);
% 
% typical post manipulations
% alpha(gca,.6);
% set(gca,'fontsize',14,'fontname','Calibri');
% ax = ancestor(gca, 'axes');
% xaxes = get(ax,'XAxis');
% set(xaxes,'fontsize',12);
% text(-.1,1.1,'A','fontsize',16,'Unit','normalized','fontname','Calibri');
% % other color-maps
% colmap = [
%     0.2000         0         0
%     0.8000    0.4000    0.1000
%     1.0000    0.2000    0.2000
%     ];


% ---------------------
if size(labels,1)>1
    ul = unique(labels(1,:),'Stable');
    mx0 = mx;
    ex0 = ex;
    labels0 = labels;
    for i=1:length(ul)
        ii = strcmp(labels(1,:),ul{i});
        m(:,i) = mx(ii);
        e(:,i) = ex(ii);        
    end
    mx = m;
    ex = e;
    glbl = unique(labels(2,:),'Stable');
    labels = ul;
end


if nargin<5, barwidth = []; end
if nargin<6, ylbl = ''; end

[K,N] = size(mx);

if size(ex,1)==K
    el    = mx-ex;
    eh    = mx+ex;
elseif size(ex,1)==(2*K)
    el    = ex(1:K,:);
    eh    = ex(K+(1:K),:);
else
    error('!');
end

kk = 0:(1/K):1; kk(end)=[];
dx = 2;

if nargin<4
if K==1
    colmap = [.5 .5 .5];
else
    colmap = repmat( (0:(1/(K-1)):1)',1,3);
end
end

basevalue = 0;
if isempty(barwidth)
    wb = 1/(K+5);
else
    wb = barwidth;
end
a     = nan(1,N);
% figure;
for i=1:N
    ax = -median(kk) + kk + +dx*(i-1);
    axs(i,:) = ax; %#ok<AGROW>
    a(i) = median(ax);
    for k=1:K
        hb(k,i) = bar(ax(k),mx(k,i),wb,'FaceColor',colmap(k,:),'EdgeColor','k','linewidth',1,'basevalue',basevalue);
        hold on;
    end
    for k=1:K        
        plot([ax(k);ax(k)],[el(k,i);eh(k,i)],'-','color','k','linewidth',2);
    end
    ylabel(sprintf('%s',ylbl));
end

if N>1
    % labels are determined by columns
    xt = a;
elseif N==1
    % labels are determined by rows
    xt = ax;
end

set(gca,'xtick',xt);
if ~isempty(labels)
    set(gca,'xticklabel',labels);
end

if K>1
dd = axs(1,2)-axs(1,1);
xlims = [axs(1,1)-dd axs(end,end)+dd];
set(gca,'xlim',xlims);
end

% set axes propertis
% set(gca,'box','off');
set(gca,'ticklength', [0 0]);


end