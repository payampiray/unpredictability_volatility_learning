function h = sim_anxiety_trait(nr,nc,subplots)

fname = fullfile(def('pipedir'),sprintf('%s.mat',mfilename));
do_sim = ~exist(fname,'file');
% do_sim = 1;

if do_sim
    [o,x,tvolatile,tstable] = timeseries;
    config = struct('v0',0.001,'u0',0.001,'nparticles',100,...
                    'tvolatile',tvolatile,'tstable',tstable,'state',x,'n',30,...
                    'rng_id',0,'nsim',100);
    rng(config.rng_id);
    nsim = config.nsim;
    n = config.n;    

    mean_triat = [.5 1 3 4 5];        
    rmin = 7^-1;
    rr = 2*(mean_triat.^-1 -rmin);
        
    lrs = cell(1,nsim);
    lambda_vs = cell(1,nsim);
    lambda_us = cell(1,nsim);
    for j=1:nsim
        lambda_v = .2*rand(1,n);
        
        K = length(rr);
        ratio = [];
        for k=1:K
            ratio = [ratio rmin+rr(k)*rand(1,6)]; %#ok<AGROW>
        end        
        lambda_u = ratio.*lambda_v;        
        config.lambda_v = lambda_v;
        config.lambda_u = lambda_u;


        N = length(o);
        outcome = nan(N,n);
        vol = nan(N,n);
        unp = nan(N,n);
        lr = nan(N,n);
        val = nan(N,n);    

        for i=1:n            
            [outcome(:,i)] = timeseries;
            conf = config;
            conf.lambda_v = lambda_v(i);
            conf.lambda_u = lambda_u(i);
            [vol(:,i),unp(:,i),lr(:,i),val(:,i)] = model_pf(outcome(:,i),conf);            
        end
        lrs{j} = lr;
        lambda_vs{j} = lambda_v;
        lambda_us{j} = lambda_u;
    
    end
    sim = struct('config',config,...
                 'lambda_vs',{lambda_vs},'lambda_us',{lambda_us},'lrs',{lrs}); %#ok<NASGU>
    save(fname,'sim');    
end

sim = load(fname); sim = sim.sim;

[trait, relative_llr] = data_anxiety_trait;
%--------------------------------------------------------------------------

t = [sim.config.tstable sim.config.tvolatile];
i = 1;
lambda = sim.lambda_vs{i}'./sim.lambda_us{i}';
lr1 = sim.lrs{i};
av = log(lr1(t(:,2),:));
as = log(lr1(t(:,1),:));
dlr = mean(av-as,1)';

nsim = length(sim.lrs);
c = nan(nsim,1);

for i=1:nsim
    lp = sim.lambda_vs{i}'./sim.lambda_us{i}';
    av = log(sim.lrs{i}(t(:,2),:));
    as = log(sim.lrs{i}(t(:,1),:));
    lrv = mean(av,1)';
    lrs = mean(as,1)';
    lr = lrv-lrs;        
    
    c(i) = corr(lp,lr,'type','spearman');  
end
mc = median(c);
ec = serr(c);
%--------------------------------------------------------------------------
if nargin<1
    close all;    
    nr = 1;
    nc = 2;
    subplots = 1:2;
    fsiz = [.3 .3 .5 .37];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

xstr = {def('lr'), def('vol'), def('unp')};

alf = def('alf');
col = def('col_br');
fsy = def('fsy');
fs = def('fs');


ylbl = sprintf('Relative log learning rate\n (volatile block - stable block)');
xlbls = {'Trait anxiety',sprintf('Model trait anxiety\n(Relative %s to %s update rate)',lower(xstr{2}),lower(xstr{3}))};
ylbl_inst = sprintf('Correlation');
ttls = {'Data','Model'};

x = {trait,lambda};
y = {relative_llr,dlr};
xl = {[15 75],[0 7]};

for i=1:2
    h(i) = subplot(nr,nc,subplots(i));
    h1 = scatter(x{i},y{i},[],col(1,:),'filled','marker','o');    
    xlim(xl{i});
    hl = lsline;
    set(hl,'color',col(1,:));
%     scatterplot(x{i},y{i},conf);    

    ylabel(ylbl,'fontsize',fsy);
    xlabel(xlbls{i},'fontsize',fsy);
    title(ttls{i},'fontsize',fsy);
end
ylim([-.5 1]);


% create smaller axes in top right, and plot on it
% axes(h(2))
h_in = axes('Position',[.8 .68 .1 .2]);
errorbarKxN(mc,ec,{''},col(2,:),.2);    
set(h_in,'box','on','LineWidth',1);

ylim([-0.5 0]);
set(gca,'fontsize',fs);
alpha(alf);
xaxes = get(gca,'XAxis');
set(xaxes,'fontsize',fsy);    
title(ylbl_inst,'fontsize',fsy,'fontweight','normal');

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
