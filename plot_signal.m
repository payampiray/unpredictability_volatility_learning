function [h, hp] = plot_signal(nr,nc,sub_plts,x,e,mstr,all_title,N,yll,abc,col,do_xlbl)

if ~iscell(sub_plts)
    sub_plts = num2cell(sub_plts);
end

if nargin<6, mstr = []; end
if nargin<7, all_title = []; end
if nargin<8, N = nan; end
if nargin<9, yll = []; end
if nargin<10, abc = ''; end
if nargin<11, col = ''; end
if nargin<12, do_xlbl = 1; end

fs = def('fs');
fn = def('fn');
fsy = def('fsy');
fsA = def('fsA');
xsA = def('xsA');
ysA = def('ysA');

linestyle = {'-';'-'};
if isempty(col)
    col = def('col');
    if size(x{1},2)>2
        col = [col;col];
        linestyle = [linestyle; '--'; '--'];
    end
end
if isempty(abc)
    abc = [];
end

if ~iscell(yll)
    yll = repmat({yll},size(x));
end

for j=1:length(x)
    hh = subplot(nr,nc,sub_plts{j});
    h(j) = hh(1);
    set(gca,'fontsize',fs,'fontname',fn);
    set_yl = 0;
    if isempty(yll{j})
        if ~isempty(e{j})
            maxy = max(max(x{j}+e{j}));
            miny = min(min(x{j}-e{j}));
        else
            maxy = max(max(x{j}));
            miny = min(min(x{j}));
        end
        
        bias = 0;
        if (maxy - miny)==0
            bias = .01;
        end
        yl = [miny maxy]+[-bias bias];
    else
        yl = yll{j};
        set_yl = 1;        
    end
    if ~isnan(N)
        yl0 = 1;
        while sum(abs(yl-yl0))>0
            yl0 = yl;
            plot([N N],yl+[10^-30 -10^-30],'color',.8*[1 1 1],'linewidth',2); hold on;
            yl = get(h(j),'ylim');
        end        
    end    
    
    for i=1:size(x{j},2)
        if isempty(e{j})
            hp(j,i) = plot(x{j}(:,i),'linewidth',2,'color',col(i,:),'linestyle',linestyle{i}); hold on;
        else
            [hp(j,i),hf(j,i)] = plot_ci(x{j}(:,i),e{j}(:,i),col(i,:)); hold on;
            set(hp(i),'linewidth',2,'linestyle',linestyle{i});
        end
        set(gca,'fontsize',fs);
    end
    nt = size(x{1},1);
    xlim([0 nt]);
    
    if set_yl
        set(h(j),'ylim',yl);
    end
    
    if ~isempty(mstr{j})
        ylabel(mstr{j},'fontsize',fsy);
    end
    
    if do_xlbl
        xlabel('Trial','fontsize',fsy);  
    end
        
    if ~isempty(abc)
        text(xsA,ysA,abc(j),'fontsize',fsA,'Unit','normalized','fontname',fn);    
    end
end

if ~isempty(all_title)
    text(h(1),-.25,.5,all_title,'fontsize',fsA,'Unit','normalized','HorizontalAlignment','Center','rotation',90);
end

end

function [h, hp] = plot1(nr,nc,sub_plts,x,e,mstr,all_title,N,yll,abc,col)
if nargin<6, mstr = []; end
if nargin<7, all_title = []; end
if nargin<8, N = nan; end
if nargin<9, yll = []; end
if nargin<10, abc = ''; end
if nargin<11, col = ''; end

fs = def('fs');
fn = def('fn');
fsy = def('fsy');
fsA = def('fsA');
xsA = def('xsA');
ysA = def('ysA');

linestyle = {'-';'-'};
if isempty(col)
    col = def('col');
    if size(x{1},2)>2
        col = [col;col];
        linestyle = [linestyle; '--'; '--'];
    end
end
if isempty(abc)
    abc = [];
end

if ~iscell(yll)
    yll = repmat({yll},size(x));
end

for j=1:length(x)
    h(j) = subplot(nr,nc,sub_plts(j));
    set(gca,'fontsize',fs,'fontname',fn);
    set_yl = 0;
    if isempty(yll{j})
        if ~isempty(e{j})
            maxy = max(max(x{j}+e{j}));
            miny = min(min(x{j}-e{j}));
        else
            maxy = max(max(x{j}));
            miny = min(min(x{j}));
        end
        
        bias = 0;
        if (maxy - miny)==0
            bias = .01;
        end
        yl = [miny maxy]+[-bias bias];
    else
        yl = yll{j};
        set_yl = 1;        
    end
    if ~isnan(N)
        yl0 = 1;
        while sum(abs(yl-yl0))>0
            yl0 = yl;
            plot([N N],yl+[10^-30 -10^-30],'color',.8*[1 1 1],'linewidth',2); hold on;
            yl = get(h(j),'ylim');
        end        
    end    
    
    for i=1:size(x{j},2)
        if isempty(e{j})
            hp(j,i) = plot(x{j}(:,i),'linewidth',2,'color',col(i,:),'linestyle',linestyle{i}); hold on;
        else
            [hp(j,i),hf(j,i)] = plot_ci(x{j}(:,i),e{j}(:,i),col(i,:)); hold on;
            set(hp(i),'linewidth',2,'linestyle',linestyle{i});
        end
        set(gca,'fontsize',fs);
    end
    nt = size(x{1},1);
    xlim([0 nt]);
    
    if set_yl
        set(h(j),'ylim',yl);
    end
    
    if ~isempty(mstr{j})
        ylabel(mstr{j},'fontsize',fsy);
    end
    xlabel('Trial','fontsize',fsy);  
        
    if ~isempty(abc)
        text(xsA,ysA,abc(j),'fontsize',fsA,'Unit','normalized','fontname',fn);    
    end
end

if ~isempty(all_title)
    text(h(1),-.25,.5,all_title,'fontsize',fsA,'Unit','normalized','HorizontalAlignment','Center','rotation',90);
end

end
