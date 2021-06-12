function h = sim_anxiety_huang(nr,nc,subplots)

fname = fullfile(def('pipedir'),sprintf('%s.mat',mfilename));
do_sim = ~exist(fname,'file');

if do_sim
    config = struct('lambda_v',.2,'lambda_u',.2,'v0',.1,'u0',.1,'nparticles',100,'u0_lesioned',0.05,'beta',3,...
                    'N',270,'rng_id',0,'nsim',100);
    rng(config.rng_id);   
    N = config.N;
    nsim = config.nsim;
    
    outcome = nan(N,nsim);
    vol = nan(N,nsim);
    unp = nan(N,nsim);
    lr = nan(N,nsim);
    val = nan(N,nsim);
    
    outcomes = cell(1,2);
    vols = cell(1,2);
    unps = cell(1,2);
    lrs = cell(1,2);
    vals = cell(1,2); 
    
    lnames = {'Healthy', sprintf('%s lesion',def('unp'))};                
    groups = {'Control','Anxious'};
    for j=1:2        
        for i=1:nsim
            [outcome(:,i)] = timeseries;         
            [vol(:,i),unp(:,i),lr(:,i),val(:,i)] = model_pf(outcome(:,i),config,lnames{j});
            
        end                
        outcomes{j} = outcome;
        vols{j} = vol;
        unps{j} = unp;
        lrs{j} = lr;
        vals{j} = val;
    end
    sim = struct('config',config,'outcomes',{outcomes},'specs',{groups},...
                 'vols',{vols},'unps',{unps},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>
    save(fname,'sim');    
end
sim = load(fname); sim = sim.sim;

[m_haung, e_haung] = data_anxiety_huang;
%--------------------------------------------------------------------------
nsim = sim.config.nsim;
ncond = length(sim.lrs);

rng(sim.config.rng_id);
beta  = sim.config.beta;

ws = nan(nsim,ncond);
ls = nan(nsim,ncond);
for j=1:ncond
    for i=1:nsim
        o = sim.outcomes{j}(:,i);
        val = sim.vals{j}(:,i);
        [ws(i,j),ls(i,j)] = WSLS(o,val,beta);      
    end
end

mlr = nan(1,ncond);
elr = nan(1,ncond);

for j=1:ncond    
    lrs = mean(sim.lrs{j},1)';            
    mlr(j) = mean(lrs);
    elr(j) = serr(lrs);        
end


mws = mean(ws);
ews = serr(ws);
mls = mean(ls);
els = serr(ls);

mw = [mws; mls];
ew = [ews; els];

mhw = [m_haung(1:2); m_haung(3:4)];
ehw = [e_haung(1:2); e_haung(3:4)];

%--------------------------------------------------------------------------
if nargin<1
    close all;
    nr = 1;
    nc = 2;
    subplots = 1:2;
    fsiz = [.3 .3 .35 .2];      
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

labels = sim.specs;
col = def('col_br');
fsy = def('fsy');
alf = def('alf');
fs = def('fs');

xstr = {def('lr'), def('vol'), def('unp')};

% h = plot_bar(nr,nc,3:4,{mlr,mvol},{elr,evol},labels,xstr(1:2));
% set(h(1),'ylim',[.5 1]);
% [hx, hp] = plot_signal(nr,nc,subplots(3:4),{ lr,vol},{e_lr,e_vol},xstr(1:2),'',nan,[],'',col);

h(1) = plot_bar(nr,nc,subplots(1),{mhw},{ehw},{'Win-stay','Lose-shift'},{'Stay/Shift probability'},'',col);
legend(h(1),labels,'fontsize',fsy,'location','northeast','box','off');
title('Data');

h(2) = plot_bar(nr,nc,subplots(2),{mw},{ew},{'Win-stay','Lose-shift'},{'Stay/Shift probability'},'',col);
% legend(h(2),labels,'fontsize',fsy,'location','northeast','box','off');
title('Model');

% h(3) = plot_bar(nr,nc,subplots(3),{mlr},{elr},labels,xstr(1),'',col);
% ylim([0 1]);


% create smaller axes in top right, and plot on it
% axes(h(2))
wp = .3;
hp = .3;
inpos = get(h(2),'position');
inpos(1) = inpos(1)+(1-1.23*wp)*inpos(3);
inpos(2) = inpos(2)+(1-1.4*hp)*inpos(4);
inpos(3) = wp*inpos(3);
inpos(4) = hp*inpos(4);

h_in = axes('Position',inpos);
errorbarKxN(mlr',elr',{''},col,.1);    
set(h_in,'box','on','LineWidth',1,'fontsize',fs);
set(h_in,'ytick',0:.4:8);
ylim([0 1]);
alpha(alf);
title(xstr{1},'fontsize',fsy,'fontweight','normal');

end

function [y1,x1]=timeseries
N = 270;
omega = 0.01;
p1 = [1 3 9 9 1 3 9 3 1 1]/13;

x1 = nan(N,1);
ilast = 0;
for i=1:10
    n0 = 30+1*randn;
    n0 = round(n0);
    ii = ilast + (1:n0);
    n  = round((1-p1(i))*n0);
    j = randperm(n0);
    x1(ii) = ones(n0,1);
    x1(ii(j(1:n)) ) = 0;
    ilast = ii(end);
end
x1 = x1(1:N);

y1 = x1 + sqrt(omega)*randn(N,1);
end

function [ws, ls, wsls] = WSLS(outcome,val,beta)

dv = 2*([0; val(1:end-1)]-.5);
p = 1./(1+exp(-beta*dv));


ps = [1 1];
choice = nan(size(p));
for t=1:length(p)
    p0 = [p(t) 1-p(t)];
    pii = p0.*ps;    
    pi = pii/sum(pii);    
    c =  binornd(1,pi(1));
    c = 2-c;    
    choice(t,1) = c;
end

r = outcome.*(choice==1)-outcome.*(choice==2);

% rpre = [nan;r(1:end-1)];
% dc = [0; diff(choice)];
% WS = (dc==0).*(rpre>.5);
% LS = (dc~=0).*(rpre<.5);

rpre = r(1:end-1);
stay = choice(2:end) == choice(1:end-1);
x = stay.*(rpre>0) + (1-stay).*(rpre<0);
wsls = mean(x);

WS = nan(size(stay));
WS(rpre>0) = stay(rpre>0);
LS = nan(size(stay));
LS(rpre<0) = 1-stay(rpre<0);

ws = nanmean(WS);
ls = nanmean(LS);
    
end