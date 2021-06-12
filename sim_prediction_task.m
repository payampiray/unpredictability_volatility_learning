function h = sim_prediction_task(nr,nc,subplots)

fname = fullfile(def('pipedir'),sprintf('%s.mat',mfilename));
do_sim = ~exist(fname,'file');

if do_sim

    config = struct('lambda_v',.4,'lambda_u',.2,'v0',5,'u0',5,'nparticles',100,...
                    'noise',[1 9],'N',300,'specs',{{'Low','High'}},...
                    'rng_id',0,'nsim',100);
    rng(config.rng_id);      
    nsim = config.nsim;
    N = config.N;
    noise = config.noise;
    
    [points] = change_points;
    state = nan(N,nsim);    
    for i=1:nsim
        [state(:,i)] = timeseries(points);
    end    
    
    ncond = length(noise);
    outcomes = cell(1,ncond);
    vols = cell(1,ncond);
    unps = cell(1,ncond);
    lrs = cell(1,ncond);
    vals = cell(1,ncond);    
    for j=1:ncond 
        outcome = nan(N,nsim);
        vol = nan(N,nsim);
        unp = nan(N,nsim);
        lr  = nan(N,nsim);
        val = nan(N,nsim);        
        for i=1:nsim            
            outcome(:,i) = state(:,i) + sqrt(noise(j))*randn(N,1);
            [vol(:,i),unp(:,i),lr(:,i),val(:,i)]=model_pf(outcome(:,i),config);
        end
        outcomes{j} = outcome;
        vols{j} = vol;
        unps{j} = unp;
        lrs{j} = lr;
        vals{j} = val;
    end
    sim = struct('config',config,'specs',{specs},...
                 'vols',{vols},'unps',{unps},'lrs',{lrs},'vals',{vals},'state',state,'outcomes',{outcomes}); %#ok<NASGU>
    save(fname,'sim');    
end
sim = load(fname); sim = sim.sim;

%--------------------------------------------------------------------------

noise = sim.config.noise;
ncond = length(noise);
nsim = size(sim.lrs{1},2);

N = sim.config.N;
state = sim.state;
tpost = 20;
tpre = 4;
t_origin_cp = (-tpre:tpost);

mlr_change_points = nan(ncond,length(t_origin_cp));
elr_change_points = nan(ncond,length(t_origin_cp));

bins_err = [0:0.05:4 inf];
lr_err_bins = nan(nsim,length(bins_err)-1);
mlr_err_bins = nan(ncond,length(bins_err)-1);
elr_err_bins = nan(ncond,length(bins_err)-1);
a = nan(nsim,ncond);
outcome = nan(N,ncond);

for j=1:ncond      
    lrs = sim.lrs{j};    
    lr_chane_points = nan(nsim,length(t_origin_cp));
    for n=1:nsim        
        points = find([0; diff(state(:,n))]~=0 );
        points(end) = [];
        lr_cp = zeros(0,length(t_origin_cp));
        for i=1:length(points)
            tcp = points(i)+t_origin_cp;        
            lr_cp(i,:) = lrs(tcp,n)';        
        end    
        lr_chane_points(n,:) = median(lr_cp,1);
        
        
        o = sim.outcomes{j}(:,n);        
        val = sim.vals{j}(:,n);
        
        val = [0; val(1:end-1)];
        delta = abs(o - val);
        
        z_del = delta/sqrt(noise(j));
        for i=2:size(bins_err,2)
            tt = (z_del>bins_err(i-1)) & (z_del<bins_err(i));
            lr_err_bins(n,i-1) = nanmedian(lrs(tt,n));
        end        
    end       
        
    mlr_change_points(j,:) = mean(lr_chane_points,1);
    elr_change_points(j,:) = serr(lr_chane_points,1);
    
    mlr_err_bins(j,:) = nanmean(lr_err_bins);
    elr_err_bins(j,:) = nanserr(lr_err_bins);
    
    st = sim.state(:,1);
    outcome(:,j) = sim.outcomes{j}(:,1);        
    a(:,j) = median(sim.lrs{j},1);    
end

mlr = mean(a);
elr = serr(a);

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    nr = 2;
    nc = 2;
    subplots = 1:4;
    fsiz = [0 0 .45 .45];    
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

levels = {'Small','Large'};
lgtitle = sprintf('Noise level');

xstr = {def('lr'), def('vol'), def('unp')};

fsy = def('fsy');
col = def('col');

yll = [-10 32];

[hx, hp] = plot_signal(nr,nc,subplots(1),{ outcome},{zeros(N,2)},{'Outcomes'},'',nan,yll,'',col);
lg = legend(hp(1,:),levels,'fontsize',fsy,'location','northeast','box','off','AutoUpdate','off','orientation','horizontal');
title(lg,lgtitle);
hold on;
plot(hx(1), outcome(:,1),'color',col(1,:),'linewidth',2);
plot(hx(1), st,'k','linewidth',2);

h = hx;

hx = plot_bar(nr,nc,subplots(2),{mlr},{elr},levels,xstr);
ylim([0 .56]);
xlabel(lgtitle,'fontsize',fsy);

h(2) = hx;

%--------------------------------------------------------------------------
levels = {'Small','Large'};
lgtitle = sprintf('Noise level');

xstr = {def('lr'), def('vol'), def('unp')};


fsy = def('fsy');
col = def('col');

i = 3;
h(i) = subplot(nr,nc,subplots(i));
for j=1:2    
    xm = mlr_change_points(j,:)';
    xe = elr_change_points(j,:)';
    plot_ci(xm,xe,col(j,:)); hold on;
    xlim([1 length(xm)]);    
end
xtl = get(gca,'xtick');
dxtl = xtl(2)-xtl(1);
xl = t_origin_cp(dxtl:dxtl:end);
set(gca,'xticklabel',xl);   
ylabel(xstr{1},'fontsize',fsy);
xlabel('Trials after change point','fontsize',fsy);  
% text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn);    

i = 4;
h(i) = subplot(nr,nc,subplots(i));
xt = 0:10:90;
for j=1:2    
    xm = mlr_err_bins(j,:)';
    xe = elr_err_bins(j,:)';
    [hp(j)]=plot_ci(xm,xe,col(j,:)); hold on;
    set(gca,'xtick',xt);
    xlim([0 length(xm)]);    
end
xtl = get(gca,'xtick');
dxtl = xtl(2)-xtl(1);

xl = bins_err(dxtl:dxtl:end);
if xtl(1)==0
    xl = [0 xl];
end
set(gca,'xticklabel',xl);   

ylabel(xstr{1},'fontsize',fsy);
xlabel('Error magnitude','fontsize',fsy);  
% text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn);    
lg = legend(hp,levels,'fontsize',fsy,'location','northwest','box','off','AutoUpdate','off');
title(lg,lgtitle);
end

function se = nanserr(x,dim)
if nargin<2, dim=1; end
s = nanstd(x,[],dim);
n = sum(~isnan(x),dim);
se = s./sqrt(n);
end

function [points] = change_points
N = 300;
rate = 0.05;
tmin = 5;

ok  = 1;
lst = 0;
points = [];
while ok
    r = tmin + ceil(exprnd(rate^-1));
    points = [points lst+r]; %#ok<AGROW>
    lst = points(end);
    ok = lst<N;
end
points(end) = [];
end

function x = timeseries(points)
N = 300;
range = 10;

x = nan(N,1);
x(1) = range*rand;
k = 2;
for t=2:N        
    jump = any(t==points);
    if jump
        x(t) = range*rand;
        k = k+1;
    else
        x(t) = x(t-1);
    end
end
end
