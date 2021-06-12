function h = sim_anxiety(nr,nc,subplots)
% anxious show insensitivity to volatility manipulation

fname = fullfile(def('pipedir'),sprintf('%s.mat',mfilename));
do_sim = ~exist(fname,'file');

if do_sim
    [o,x,tvolatile,tstable] = timeseries;
    config = struct('lambda_v',.2,'lambda_u',.2,'v0',.1,'u0',.1,'nparticles',100,'u0_lesioned',0.001,...
                    'tvolatile',tvolatile,'tstable',tstable,'state',x,...
                    'rng_id',0,'nsim',100);
    rng(config.rng_id); 
    nsim = config.nsim;
    
    N = length(o);
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
    glabels = {'Control','Anxious'};
    lnames = {'Healthy', sprintf('%s lesion',def('unp'))};            
    for j=1:2
        for i=1:nsim            
            [o] = timeseries;
            [vol(:,i),unp(:,i),lr(:,i),val(:,i)] = model_pf(o,config,lnames{j});            
        end        
        if j==1
            o_example = o;
        end
        v_example(:,j) = val(:,1);
        
        outcomes{j} = outcome;
        vols{j} = vol;
        unps{j} = unp;
        lrs{j} = lr;
        vals{j} = val;
    end
    specs(1,:) = glabels;
    
    sim = struct('config',config,'specs',{specs},...
                 'o_example',o_example,'v_example',v_example,...
                 'vols',{vols},'unps',{unps},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>
    save(fname,'sim');    
end

sim = load(fname); sim = sim.sim;

%--------------------------------------------------------------------------
val = sim.v_example;
state = sim.config.state;
N = length(sim.config.state);

glabels = sim.specs;
t = [sim.config.tstable sim.config.tvolatile];
strcols = {'Stable','Volatile'};

vol = nan(N,2);
unp = nan(N,2);
lr = nan(N,2);
evol = nan(N,2);
eunp = nan(N,2);
elr = nan(N,2);

vols = cell(N,2);
unps = cell(N,2);
lrs = cell(N,4);
for j=1:2
    vol(:,j) = mean(sim.vols{j},2);
    unp(:,j) = mean(sim.unps{j},2);
    lr(:,j) = mean(sim.lrs{j},2);
    evol(:,j) = serr(sim.vols{j},2);
    eunp(:,j) = serr(sim.unps{j},2);
    elr(:,j) = serr(sim.lrs{j},2);    
    
    for k=1:2
        l = (j-1)*2 + k;

        vols{l} = sim.vols{j}(t(:,k),:);
        unps{l} = sim.unps{j}(t(:,k),:);
        lrs{l} = sim.lrs{j}(t(:,k),:);

        vals{l} =  sim.vals{j}(t(:,k),:);

        specs{1,l} = glabels{j};
        specs{2,l} = strcols{k};            
        specs{3,l} = sprintf('%s-%s',glabels{j},strcols{k});            
    end
end

ncond = 4;
mvol = nan(1,ncond);
munp = nan(1,ncond);
mlr =  nan(1,ncond);
evol = nan(1,ncond);
eunp = nan(1,ncond);
elr =  nan(1,ncond);
for j=1:ncond
    v = mean(vols{j},2);
    u = mean(unps{j},2);
    a = mean(lrs{j},2);    
    
    mvol(j) = mean(v);
    munp(j) = mean(u);
    mlr(j) = mean(a);
    evol(j) = 0*serr(v);
    eunp(j) = 0*serr(u);
    elr(j) = serr(a);
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    nr = 2;
    nc = 2;
    subplots = 1:4;
    fsiz = [.3 .3 .35 .45];    

    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

labels = specs(1:2,:);
clabels = specs(2,[1 2]);
xstr = {def('lr'), def('vol'), def('unp')};

alf = .3;
col = def('col_br');
fsy = def('fsy');


N = size(val,1);
Nline = nan;
ii = [2 1];
[hx, hp] = plot_signal(nr,nc,subplots(1),{ val(:,ii)},{zeros(N,2)},{'Estimated reward'},'',Nline,[],'',col(ii,:));
hold on;
plot(hx, state,'color',.6*ones(1,3),'linewidth',2);

h(1) = hx;

h(2) = plot_bar(nr,nc,subplots(2),{mlr},{elr},labels,xstr(1));
lg = legend(h(2),clabels,'fontsize',fsy,'location','north','box','off');

ii = [1 2];
[hx, hp] = plot_signal(nr,nc,subplots(3:4),{ vol(:,ii), unp(:,ii)},{evol(:,ii), eunp(:,ii)},xstr(2:3),'',nan,[],'',col);
lg = legend(hp(1,:),glabels,'fontsize',fsy,'location','northwest','box','off','autoupdate','off');

h(3:4) = hx;

%--------------------------------------------------------------------------
defcol = def('col');
tstable  = find(sim.config.tstable');
tvolatile  = [tstable(end)+.5 find(sim.config.tvolatile')];
epsil = .001;

yls = {[-.09 1.09],nan,[0 .11],[0 .11]};

for i=[1 3 4]
    hs(i) = subplot(nr,nc,subplots(i));
    axes(hs(i)); %#ok<LAXES>
    hold on;
    yl = get(hs(i),'ylim')+[epsil -epsil];
    yl = yls{i};
    
    x = tstable;
    x2 = [x, fliplr(x)];
    inBetween = [yl(1)*ones(1,length(x)), yl(2)*ones(1,length(x))];
    fill(x2, inBetween, defcol(1,:), 'FaceAlpha', alf, 'EdgeColor', defcol(1,:),'EdgeAlpha', alf); hold on;     
    x = tvolatile;
    x2 = [x, fliplr(x)];
    inBetween = [yl(1)*ones(1,length(x)), yl(2)*ones(1,length(x))];
    fill(x2, inBetween, defcol(2,:), 'FaceAlpha', alf, 'EdgeColor', defcol(2,:),'EdgeAlpha', alf); hold on; 
    ylim(yl);    
end

end

function [y,x,tvolatile,tstable]=timeseries
n = 20;
m = [.75 .75*ones(1,5) .2 .8 .2 .8 .2];
omega = .01;

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
