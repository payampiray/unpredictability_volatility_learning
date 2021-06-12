function h = sim_2x2(nr,nc,subplots,fig_no)

factors.unp = {'Small','Small','Large','large'};
factors.vol = {'Small','Large','Large','Small'};
fname = fullfile(def('pipedir'),sprintf('%s.mat',mfilename));
do_sim = ~exist(fname,'file');

if do_sim
    config = struct('x0_unc',100,'lambda_v',.1,'lambda_u',.1,'v0',1,'u0',2,'nparticles',100,...
                    'true_vol',[.5 1.5 .5 1.5],'true_unp',[1 1 3 3],'factors',factors,'N',200,...
                    'rng_id',0,'nsim',100);
                
    rng(config.rng_id);
    true_vol = config.true_vol;
    true_unp = config.true_unp;
    N = config.N;
    nsim = config.nsim;    
    
    outcomes = cell(1,4);
    vols = cell(1,4);
    unps = cell(1,4);
    lrs = cell(1,4);
    vals = cell(1,4);
    
    specs = cell(3,4);
    for j=1:4
        outcome = nan(N,nsim);        
        vol = nan(N,nsim);
        unp = nan(N,nsim);
        lr  = nan(N,nsim);
        val = nan(N,nsim);    
        for i=1:nsim
            [outcome(:,i)] = timeseries(config,true_vol(j),true_unp(j));        
            [vol(:,i),unp(:,i),lr(:,i),val(:,i)]=model_pf(outcome(:,i),config);
        end
        outcomes{j} = outcome;        
        specs{1,j} = sprintf('%s %s',factors.unp{j},lower(def('unp')));
        specs{2,j} = sprintf('%s %s',factors.vol{j},lower(def('vol')));
        specs{3,j} = sprintf('true_vol=%0.2f, true_unp=%0.2f',true_vol(j),true_unp(j));        
        
        vols{j} = vol;
        unps{j} = unp;
        lrs{j} = lr;
        vals{j} = val;
    end
    config.specs = specs;    
    sim = struct('config',config,'specs',{config.specs},'outcomes',{outcomes},...
                 'vols',{vols},'unps',{unps},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>
    save(fname,'sim');
end
sim = load(fname); sim = sim.sim;

%--------------------------------------------------------------------------

N = sim.config.N;
nsim = sim.config.nsim;

val = nan(N,4);
vol = nan(N,4);
unp = nan(N,4);
lr =  nan(N,4);
e_val = nan(N,4);
e_vol = nan(N,4);
e_unp = nan(N,4);
e_lr =  nan(N,4);

tend = (.9*N+1):N;
a = nan(nsim,4);
v = nan(nsim,4);
u = nan(nsim,4);

for j=1:4    
    val(:,j) = mean(sim.vals{j},2);
    vol(:,j) = mean(sim.vols{j},2);
    unp(:,j) = mean(sim.unps{j},2);
    lr(:,j) = mean(sim.lrs{j},2);
    
    e_val(:,j) = serr(sim.vals{j},2);
    e_vol(:,j) = serr(sim.vols{j},2);
    e_unp(:,j) = serr(sim.unps{j},2);
    e_lr(:,j) = serr(sim.lrs{j},2);
    
    a(:,j) = mean(sim.lrs{j}(tend,:),1);
    v(:,j) = mean(sim.vols{j}(tend,:),1);
    u(:,j) = mean(sim.unps{j}(tend,:),1);
    
end

mvol = mean(v);
munp = mean(u);
mlr = mean(a);
evol = serr(v);
eunp = serr(u);
elr = serr(a);

ii1 = [1 3];
ii2 = [2 4];
mlr = [mlr(ii1)' mlr(ii2)'];
mvol = [mvol(ii1)' mvol(ii2)'];
munp = [munp(ii1)' munp(ii2)'];

elr = [elr(ii1)' elr(ii2)'];
evol = [evol(ii1)' evol(ii2)'];
eunp = [eunp(ii1)' eunp(ii2)'];

%--------------------------------------------------------------------------
if nargin<1
    close all;
    nr = 3;
    nc = 2;     
    subplots = 1:6;    
    fsiz = [0 0 .6 .4];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
    fig_no = 1;    
end
colstrs = {'Small','Large'};
glbl = {sprintf('Small true %s',lower(def('vol'))),sprintf('Large true %s',lower(def('vol')))};
levels = {'Small','Large'};
lgtitle = sprintf('True %s',lower(def('unp')));
xltitle = sprintf('True %s',lower(def('vol')));

xstr = {def('lr'), def('vol'), def('unp')};

ylb = {[0 .75],[0 2],[0 3]};
yls = {[0.2 .81],[0 1.8],[0 3.5]};

fsy = def('fsy');
abc = def('abc');

if fig_no==1
    ii = ii1;
    [hx, hp] = plot_signal(nr,nc,subplots(1:2:6),{lr(:,ii),vol(:,ii),unp(:,ii)},...
                                               {e_lr(:,ii),e_vol(:,ii),e_unp(:,ii)},xstr,'',nan, yls, '', [], 0);
    title(hx(1),glbl{1},'fontsize',fsy);
    xlabel(hx(3),'Trial','fontsize',fsy);
    lg = legend(hp(2,:),levels,'fontsize',fsy,'box','off','orientation','horizontal');    
    title(lg,lgtitle);

    h(1:3) = hx;

    ii = ii2;
    [hx] = plot_signal(nr,nc,subplots(2:2:6),{lr(:,ii),vol(:,ii),unp(:,ii)},...
                                               {e_lr(:,ii),e_vol(:,ii),e_unp(:,ii)},xstr,'',nan, yls, '', [], 0);
    title(hx(1),glbl{2},'fontsize',fsy);
    xlabel(hx(3),'Trial','fontsize',fsy);
    h(4:6) = hx;
end

%--------------------------------------------------------------------------
if fig_no == 2
    nr = 1;
    nc = 3;     
    subplots = 1:3;    
    fsiz = [0 0 .6 .25];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

    h(1:3) = plot_bar(nr,nc,subplots(1:3),{mlr',mvol',munp'},{elr',evol',eunp'},colstrs,xstr,abc,[],ylb);
    loc = {'northwest','northwest','north'};
    for i=1:length(h)
        xlabel(h(i),xltitle,'fontsize',fsy);
        if i==2
        lg = legend(h(i),levels,'fontsize',fsy,'location',loc{i},'box','off');
        title(lg,lgtitle);    
        end
    end
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
