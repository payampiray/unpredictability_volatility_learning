function h = sim_feasibility(nr,nc,subplots,fig_no)

factors.vol = {'Small','Large','Small','Large'};
factors.unp = {'Small','Small','Large','large'};
fname = fullfile(def('pipedir'),sprintf('%s.mat',mfilename));
do_sim = ~exist(fname,'file');
if do_sim
    config = struct('true_vol',[.5 1.5 .5 1.5],'true_unp',[1 1 3 3],'factors',factors,'N',200,...
                    'rng_id',10,'nsim',100);
    rng(config.rng_id);
    
    true_vol = config.true_vol;
    true_unp = config.true_unp;
    N = config.N;
    nsim = config.nsim;
    
    outcome = nan(N,4);
    state = nan(N,4);
    specs = cell(3,4);
    for j=1:4
        if j<=2
            [outcome(:,j),state(:,j)] = timeseries(config,true_vol(j),true_unp(j));
        else
            state(:,j) = state(:,j-2);
            outcome(:,j) = state(:,j) + sqrt(true_unp(j))*randn(N,1);
        end        
        specs{1,j} = sprintf('%s %s',factors.unp{j},lower(def('unp')));
        specs{2,j} = sprintf('%s %s',factors.vol{j},lower(def('vol')));
        specs{3,j} = sprintf('$\\sigma=%0.2f$, $\\omega=%0.2f$',true_vol(j),true_unp(j));        
    end
    ac1 = nan(nsim,4);
    var1 = nan(nsim,4);
    for j=1:4        
        for i=1:nsim
            [o] = timeseries(config,true_vol(j),true_unp(j));
            [a, lags] = autocorr(o,'NumLags',1); %#ok<ASGLU>
            ac1(i,j) = a(2);
            var1(i,j) = std(o);            
        end
    end    
    sim = struct('config',config,'specs',{specs},...
                 'outcome',outcome,'state',state,'var1',var1,'ac1',ac1); %#ok<NASGU>    
    save(fname,'sim');     
end
sim = load(fname); sim = sim.sim;

%--------------------------------------------------------------------------

if nargin< 1
    close all;    
    nr = 2;
    nc = 1;
    fsiz = [0 0 .4 .4];
    subplots = 1:2;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
    fig_no = 1;
end


outcome = sim.outcome;
state = sim.state;
mv = mean(sim.var1);
ev = serr(sim.var1);

mac1 = mean(sim.ac1);
eac1 = serr(sim.ac1);

mac1 = [mac1([1 3])' mac1([2 4])'];
eac1 = [eac1([1 3])' eac1([2 4])'];

mv = [mv([1 3])' mv([2 4])'];
ev = [ev([1 3])' ev([2 4])'];

colstrs = {'Small','Large'};
glbl = {'Small','Large'};
lgtitle = def('unp');
olbl = sim.specs(2,[1 2]);
ii = [1 3;2 4];


col = def('col');
fn = def('fn');
fsy = def('fsy');
fsl = fsy+2;
yl = [-20 20];

if fig_no == 1
    for i=1:2
        h(i) = subplot(nr,nc,subplots(i));
        hp = nan(1,2);    
        hx = nan(1,2);
        plot(state(:,i),'-','linewidth',2,'color','k');

        hold on;
        for j=1:2
            ij = ii(i,j);
            hp(j) = plot(outcome(:,ij),'-','color',col(j,:),'linewidth',1); hold on;
            hx(j) = plot(outcome(:,ij),'.','color',col(j,:),'markersize',10); hold on;
        end    
        set(gca,'fontname',fn,'box','off');
        title(olbl{i},'fontsize',fsl);     
        if i==1
            lg = legend(hp,glbl,'fontsize',fsy,'box','off','orientation','horizontal');    
            title(lg,lgtitle);
        end
        ylim(yl);        
        ylabel('Outcome','fontsize',fsy);
        if i==2
            xlabel('Trial','fontsize',fsy);
        end
    end
end

%--------------------------------------------------------------------------
if nargin<1
    nr = 1;
    nc = 2;
    fsiz = [0 0 .4 .25];
    subplots = 1:2;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
    fig_no = 2;
end

if fig_no == 2
    hx = plot_bar(nr,nc,subplots(1:2),{mv',mac1'},{ev',eac1'},colstrs,{'Variance','Autocorrelation'});
    lg = legend(hx(1),glbl,'fontsize',fsy,'location','northwest','box','off');
    title(lg,lgtitle);
    for i=1:2
        xlabel(hx(i),def('vol'),'fontsize',fsy); 
    end
    set(hx(2),'ylim',[.5 1]);
end

end

function [y,x]=timeseries(config,true_vol,true_unp)
N = config.N;
x = zeros(N,1);
for t=2:N
    x(t) = x(t-1) + sqrt(true_vol)*randn;
end
y = x + sqrt(true_unp)*randn(N,1);
end

