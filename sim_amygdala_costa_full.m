function h = sim_amygdala_costa_full(nr,nc,subplots)

fname = fullfile(def('pipedir'),sprintf('%s.mat',mfilename));
do_sim = ~exist(fname,'file');

if do_sim
    config = struct('lambda_v',.1,'lambda_u',.1,'v0',.5,'u0',.5,'nparticles',100,...
                    'v0_lesioned',0.01,'beta',3,'beta_lesioned',1,...
                    'omega',10^-6,'p',[.6 .7 .8 1],'N',80,...
                    'rng_id',0,'nsim',100);
    rng(config.rng_id);
    nsim = config.nsim;
    N = config.N;
    p = config.p;
    J = 4;
    
    outcomes = cell(1,2*J);
    vols = cell(1,2*J);
    unps = cell(1,2*J);
    lrs = cell(1,2*J);
    vals = cell(1,2*J);
    k = 0;
    lnames = {'Control','Lesioned'};
    lesion_model = {'Healthy', sprintf('%s lesion',def('vol'))};                
    
    for l=1:2
        for j=1:J   
            jname = sprintf('%d/%d',p(j)*100,(1-p(j))*100);
            
            outcome = nan(N,nsim);
            vol = nan(N,nsim);
            unp = nan(N,nsim);
            lr  = nan(N,nsim);
            val = nan(N,nsim);
            for i=1:nsim                
                outcome(:,i) = timeseries(config,p(j));
                [vol(:,i),unp(:,i),lr(:,i),val(:,i)]=model_pf(outcome(:,i),config,lesion_model{l});
            end
            k = k+1;
            specs(:,k) = [{sprintf('%s-%s',lnames{l},jname)} lnames(l) jname sprintf('%d/%d',p(j)*100,(1-p(j))*100) k]; %#ok<AGROW>
            
            outcomes{k} = outcome;        
            vols{k} = vol;
            unps{k} = unp;
            lrs{k} = lr;
            vals{k} = val;
        end
    end
    sim = struct('config',config,'specs',{specs},...
                 'vols',{vols},'unps',{unps},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>    
    save(fname,'sim');    
end
sim = load(fname); sim = sim.sim;

data = data_amygdala_costa(1);

%--------------------------------------------------------------------------
if nargin<1
%     close all;    
    nr = 1;
    nc = 2;            
    fsiz = [0 0 .45 .3];
    subplots = 1:2;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

specs = sim.specs;
np = length(sim.config.p);
beta = [sim.config.beta*ones(1,np) sim.config.beta_lesioned*ones(1,np)];

N = sim.config.N;
N1 = N/2;
N2 = N/2;
vol = nan(N,2);
unp = nan(N,2);
lr = nan(N,2);
e_vol = nan(N,2);
e_unp = nan(N,2);
e_lr = nan(N,2);

val = nan(N,2);
e_val = nan(N,2);
for j=1:length(sim.vols)
    
    v   = sim.vals{j};
    dv  = v;
    v   = 1./(1+exp(-beta(j)*dv));
    val(:,j) = mean(v,2);
    e_val(:,j) = serr(v,2);
    
    vol(:,j) = mean(sim.vols{j},2);
    unp(:,j) = mean(sim.unps{j},2);
    lr(:,j) = mean(sim.lrs{j},2);    
        
    e_vol(:,j) = serr(sim.vols{j},2);
    e_unp(:,j) = serr(sim.unps{j},2);
    e_lr(:,j) = serr(sim.lrs{j},2);    
end

pnames = {'60/40','70/30','80/20','100/0'};

gnames = specs(2,[1 3]);
% cnames = {'Acquisition','Reversal'};
cnames = {'Acq','Rev'};


tt1 = 1:N1;
tt2 = (N1+1):(N1+N2);
mval_acq = mean(val(tt1,:));
eval_acq = serr(val(tt1,:));

mval_rev = 1-mean(val(tt2,:));
eval_rev = serr(1-val(tt2,:));

mdat_acq  = data.x_acq;
edat_acq  = data.e_acq;

mdat_rev  = data.x_rev;
edat_rev  = data.e_rev;

%-----------

h(1) = subplot(nr,nc,subplots(1));
ylbl = 'Fraction of correct choice';
plot_fraction(mdat_acq,mdat_rev,edat_acq,edat_rev,gnames,cnames,pnames,ylbl,1);
% % text(xsA,ysA,abc(2),'fontsize',fsA,'Unit','normalized','fontname',fn);        
title('Data');

h(2) = subplot(nr,nc,subplots(2));
ylbl = 'Probability of correct choice';
plot_fraction(mval_acq,mval_rev,eval_acq,eval_rev,gnames,cnames,pnames,ylbl,1);
% % text(xsA,ysA,abc(3),'fontsize',fsA,'Unit','normalized','fontname',fn);        
title('Model');

end

function [h,he, col]=plot_fraction(mval_acq,mval_rev,eval_acq,eval_rev,gnames,cnames,pnames,ylbl,do_leg)
np = 4; %length(mval_acq)/2;

col = [0 0 0];
col(2,:) = [255 85 155]/255;

fs = def('fs');
fn = def('fn');
fsy = def('fsy');

for l=1:2
    ll = (l-1)*np+(1:np);
    x = {mval_acq(ll), mval_rev(ll)};
    e = {eval_acq(ll), eval_rev(ll)};
    linetype = {'-','--'};
    for j=1:2        
        h(l,j) = plot(x{j},linetype{j},'color',col(l,:),'linewidth',2); hold on;
        he(l,j) = errorbar(x{j},e{j},linetype{j},'color',col(l,:),'linewidth',2); hold on;
        gcnames{l,j} = sprintf('%s-%s',gnames{l},cnames{j});
    end    
    xlim([1-.75 np+.75]);    
    set(gca,'xtick',1:np,'xticklabel',pnames);       
end
set(gca,'ylim',[.5 1],'ytick',0.5:.1:1,'fontsize',fs,'fontname',fn);
ylabel(ylbl,'fontsize',fsy);
xaxes = get(ancestor(gca, 'axes'),'XAxis');
set(xaxes,'fontsize',fsy,'ticklength', [0 0]);
% lg = legend(h(:,1),gnames,'fontsize',fsy,'location','northwest');

hh = [h(1:2) h(3:4)];
if do_leg
legend(hh,gcnames,'location','northwest','fontsize',fsy,'box','off');
end
end

function [y]=timeseries(config,p)
N = config.N;
omega = config.omega;

N = N/2;
n = p*N;
x01 = [ones(1,n) zeros(1,N-n)]';
x02 = 1-x01;

j = randperm(N);
x1 = x01(j);

j = randperm(N);
x2 = x02(j);

x = [x1; x2];
x(x==0) = -1;

N = length(x);
y = x + sqrt(omega)*randn(N,1);

end
