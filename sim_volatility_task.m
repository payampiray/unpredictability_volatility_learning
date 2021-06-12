function h = sim_volatility_task(nr,nc,subplots)

fname = fullfile(def('pipedir'),sprintf('%s.mat',mfilename));
do_sim = ~exist(fname,'file');

if do_sim
    [o,x,tvolatile, tstable] = timeseries;
    config = struct('lambda_v',.2,'lambda_u',.2,'v0',.1,'u0',.1,'nparticles',100,...
                    'tvolatile',tvolatile,'tstable',tstable,...
                    'rng_id',0,'nsim',100);    
    rng(config.rng_id);   
    nsim = config.nsim;
    
    N = length(o);
    outcome = nan(N,nsim);
    vol = nan(N,nsim);
    unp = nan(N,nsim);
    lr = nan(N,nsim);
    val = nan(N,nsim);
    
    for i=1:nsim
        [outcome(:,i)] = timeseries;
        [vol(:,i),unp(:,i),lr(:,i),val(:,i)] = model_pf(outcome(:,i),config);
    end        
    t = [tstable tvolatile];
    outcomes = cell(1,2);
    vols = cell(1,2);
    unps = cell(1,2);
    lrs = cell(1,2);
    vals = cell(1,2);
    for j=1:2
        outcomes{j} = o(t(:,j),:);
        vols{j} = vol(t(:,j),:);
        unps{j} = unp(t(:,j),:);
        lrs{j} = lr(t(:,j),:);
        vals{j} = val(t(:,j),:);
    end
    columns = {'Stable','Volatile'};    
    sim = struct('config',config,'specs',{columns},...
                 'vols',{vols},'unps',{unps},'lrs',{lrs},'vals',{vals},...
                 'vol',{vol},'unp',{unp},'lr',{lr},'val',{val},'outcome',{outcome},'state',x); %#ok<NASGU>
    save(fname,'sim');    
end
sim = load(fname); sim = sim.sim;

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    nr = 2;
    nc = 2;        
    fsiz = [0 0 .5 .5];
    subplots = 1:4;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

ncond = 2;
vol = cell(1,ncond);
unp = cell(1,ncond);
lr  = cell(1,ncond);
for j=1:ncond
    vol{j} = mean((sim.vols{j}),2);
    unp{j} = mean((sim.unps{j}),2);
    lr{j} = mean((sim.lrs{j}),2);
end

mvol = nan(1,ncond);
munp = nan(1,ncond);
mlr =  nan(1,ncond);
evol = nan(1,ncond);
eunp = nan(1,ncond);
elr =  nan(1,ncond);
for j=1:ncond
    mvol(j) = mean(vol{j});
    munp(j) = mean(unp{j});
    mlr(j) = mean(lr{j});
    evol(j) = serr(vol{j});
    eunp(j) = serr(unp{j});
    elr(j) = serr(lr{j});    
end

state = sim.state;
val = mean(sim.val,2);
vol = mean(sim.vol,2);
unp = mean(sim.unp,2);
e_val = serr(sim.val,2);
e_vol = serr(sim.vol,2);
e_unp = serr(sim.unp,2);

labels = sim.specs;
xstr = {def('lr'), def('vol'), def('unp')};


alf = .3;
col = def('col');

[hs(1)] = plot_signal(nr,nc,subplots(1),{val},{e_val},{'Estimated reward'},[],[],[],[],[0 0 0]);
plot(hs(1),state,'--k','linewidth',1);
ylim([0 1]);

h(2) = plot_bar(nr,nc,subplots(2),{mlr},{elr},labels,xstr(1));

[hs(2:3)] = plot_signal(nr,nc,subplots(3:4),{vol,unp},{e_vol,e_unp},xstr(2:3),[],[],[],[],[0 0 0]);

h([1 3 4]) = hs;

tstable  = find(sim.config.tstable');
tvolatile  = [tstable(end)+.5 find(sim.config.tvolatile')];
for i=1:length(hs)
    yl = get(hs(i),'ylim')+[eps -eps];
    axes(hs(i)); %#ok<LAXES>
    hold on;

    x = tstable;
    x2 = [x, fliplr(x)];
    inBetween = [yl(1)*ones(1,length(x)), yl(2)*ones(1,length(x))];
    fill(x2, inBetween, col(1,:), 'FaceAlpha', alf, 'EdgeColor', col(1,:),'EdgeAlpha', alf); hold on;     
    x = tvolatile;
    x2 = [x, fliplr(x)];
    inBetween = [yl(1)*ones(1,length(x)), yl(2)*ones(1,length(x))];
    fill(x2, inBetween, col(2,:), 'FaceAlpha', alf, 'EdgeColor', col(2,:),'EdgeAlpha', alf); hold on; 
    ylim(yl);    
end
end

function [y,x,tvolatile,tstable]=timeseries
n = 20;
m = [.75 .75*ones(1,5) .2 .8 .2 .8 .2];
omega = 0.01;

N = length(m)*n;

x = nan(N,1);
for i=1:length(m)
    ii = (i-1)*n+ (1:n);
    x(ii) = m(i)*ones(n,1);    
end

tvolatile = zeros(N,1);
tstable = zeros(N,1);
tstable(1:120) = 1;
tvolatile(121:N) = 1;
tstable(1:n) = 0;

tvolatile = tvolatile == 1;
tstable = tstable == 1;

y = x + sqrt(omega)*randn(N,1);
end
