function h = sim_lesioned(nr,nc,subplots)

fname = fullfile(def('pipedir'),sprintf('%s.mat',mfilename));
fsim = fullfile(def('pipedir'),sprintf('sim_2x2.mat'));

do_sim = ~exist(fname,'file');

if do_sim
    sim2x2 = load(fsim); sim2x2 = sim2x2.sim;
    config = sim2x2.config;
    outcomes = sim2x2.outcomes;
    
    nsim = config.nsim;
    config.v0_lesioned = config.v0;
    config.u0_lesioned = config.u0;
    
    N = config.N;    

    lnames = {'Healthy'; sprintf('%s lesion',def('unp')); sprintf('%s lesion',def('vol'))};        
    config.lnames = lnames;
    
    vols = cell(length(lnames),4);
    unps = cell(length(lnames),4);
    lrs  = cell(length(lnames),4);
    vals = cell(length(lnames),4);   
    
    vols(1,:) = sim2x2.vols;
    unps(1,:) = sim2x2.unps;
    lrs(1,:)  = sim2x2.lrs;
    vals(1,:) = sim2x2.vals;

    for l=2:3
        rng(config.rng_id);
        
        vol = nan(N,nsim);
        unp = nan(N,nsim);
        lr  = nan(N,nsim);
        val = nan(N,nsim);            
        for j=1:4
            for i=1:nsim
                [vol(:,i),unp(:,i),lr(:,i),val(:,i)]=model_pf(outcomes{j}(:,i),config,lnames{l});
            end
            vols{l,j} = vol;
            unps{l,j} = unp;
            lrs{l,j} = lr;
            vals{l,j} = val;
        end
    end
    
    sim = struct('config',config,'specs',{config.specs},...                                
                 'vols',{vols},'unps',{unps},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>
    save(fname,'sim');
end
sim = load(fname); sim = sim.sim;
%--------------------------------------------------------------------------

N = sim.config.N;
[L,J] = size(sim.lrs);
nsim = sim.config.nsim;

ttend = (.9*N+1):N;
mlr = cell(1,L);
mvol = cell(1,L);
munp = cell(1,L);
elr = cell(1,L);
evol = cell(1,L);
eunp = cell(1,L);

ii1 = [1 3];
ii2 = [2 4];

for l=1:L
    a = nan(nsim,J);
    v = nan(nsim,J);
    u = nan(nsim,J);
    for j=1:J
        a(:,j) = mean(sim.lrs{l,j}(ttend,:),1);
        v(:,j) = mean(sim.vols{l,j}(ttend,:),1);
        u(:,j) = mean(sim.unps{l,j}(ttend,:),1);
    end
    ma = mean(a);
    mv = mean(v);
    mu = mean(u);    
    
    ea = serr(a);
    ev = serr(v);
    eu = serr(u);    

    mlr{l} = [ma(ii1)' ma(ii2)']';
    mvol{l} = [mv(ii1)' mv(ii2)']';
    munp{l} = [mu(ii1)' mu(ii2)']';

    elr{l} = [ea([1 2])' ea([3 4])']';
    evol{l} = [ev([1 2])' ev([3 4])']';
    eunp{l} = [eu([1 2])' eu([3 4])']';    
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    nr = 2;
    nc = 3;
    subplots = 1:9;
    fsiz = [0 0 .7 .55];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

colstrs = {'Small','Large'};
levels = {'Small','Large'};
lgtitle = sprintf('True %s',lower(def('unp')));
xltitle = sprintf('True %s',lower(def('vol')));

lnames = sim.config.lnames;
xstr = {def('lr'), def('vol'), def('unp')};

ylb = {[0 .75],[0 .75],[0 .75]};

fsy = def('fsy');

ii = [1 2 3];
lnames = lnames(ii);
hx = plot_bar(nr,nc,subplots(1:3),mlr(ii),elr(ii),colstrs,repmat(xstr(1),1,3),'',[],ylb);
for j=1:3
    title(hx(j),lnames{j},'fontsize',fsy);
    xlabel(hx(j),xltitle,'fontsize',fsy);
end
lg = legend(hx(1),levels,'fontsize',fsy,'location','northwest','box','off');
title(lg,lgtitle);
pos = lg.Position;
pos(1) = pos(1)*1.1;
set(lg,'position',pos);
h(1:3) = hx;

hx = plot_bar(nr,nc,subplots(5:6),{mvol{2},munp{3}},{evol{2},eunp{3}},colstrs,xstr([2 3]));
for j=1:2
    title(hx(j),lnames{j+1},'fontsize',fsy);
    xlabel(hx(j),xltitle,'fontsize',fsy);
end

h(5:6) = hx;
end

