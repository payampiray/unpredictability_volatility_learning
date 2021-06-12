function h = sim_serial_prediction(nr,nc,subplots)

fname = fullfile(def('pipedir'),sprintf('%s.mat',mfilename));
do_sim = ~exist(fname,'file');
if do_sim
    config = struct('lambda_v',[.2 .2],'lambda_u',[.2 .2],'v0',[.5 .5],'u0',[.5 .5],'nparticles',100,'v0_lesioned',.25*10^-6*[1 1],...
                    'N1',200,'N2',100,'omega',10^-6,...
                    'rng_id',0,'nsim',100);
    rng(config.rng_id);
    specs = cell(4,4);
    nsim = config.nsim;
    
    N = config.N1 + config.N2;    
    outcomes = cell(1,4);
    stimuli = cell(1,4);
    vols = cell(1,4);
    unps = cell(1,4);
    lrs = cell(1,4);
    vals = cell(1,4);
    for l= [1 2]
        lnames = {'Control','Lesioned'};
        for j= [1 2]
            jnames = {'Consistent','Shift'};
            
            outcome = nan(N,nsim);  
            stimulus = nan(N,nsim);  
            vol = cell(1,2);
            unp = cell(1,2);
            lr  = cell(1,2);
            val1 = cell(1,2);
            for i=1:nsim
                [outcome(:,i), stimulus(:,i)] = timeseries(config,jnames{j});                
                o = [stimulus(:,i), outcome(:,i)];
                [v,u,a,m]=model_pfm(o,config,l==2);                
                for k=1:2
                    vol{k}(:,i) = v(:,k);
                    unp{k}(:,i) = u(:,k);
                    lr{k}(:,i) = a(:,k);
                    val1{k}(:,i) = m(:,k);
                end
            end
            k = (l-1)*2 + j;
            specs(:,k) = [{sprintf('%s-%s',lnames{l},jnames{j})} lnames(l) jnames(j) k];

            stimuli{k} = stimulus;
            outcomes{k} = outcome;        
            vols{k} = vol;
            unps{k} = unp;
            lrs{k} = lr;
            vals{k} = val1;
        end
    end
    sim = struct('config',config,'specs',{specs},'stimuli',{stimuli},'outcomes',{outcomes},...
                 'vols',{vols},'unps',{unps},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>    
    save(fname,'sim');    
end
sim = load(fname); sim = sim.sim;

[m_emp, e_emp, ph1_emp] = data_serial_prediction;
%--------------------------------------------------------------------------


specs = sim.specs;
N1 = sim.config.N1;
N2 = sim.config.N2;
nsim = sim.config.nsim;

t12 = 1:(N1+N2);
ttend = (N1+N2);

N  = length(t12);
vol = nan(N,4);
unp = nan(N,4);
lr = nan(N,4);
e_vol = nan(N,4);
e_unp = nan(N,4);
e_lr = nan(N,4);

val1 = nan(nsim,4);

dim = 1:2;
a = cell(1,2);
for k=1:length(dim)
    dimk = dim(k);
    for j=1:4
        if dimk==1
            vol(:,j) = mean(sim.vols{j}{dimk}(t12,:),2);
            unp(:,j) = mean(sim.unps{j}{dimk}(t12,:),2);
            e_vol(:,j) = serr(sim.vols{j}{dimk}(t12,:),2);
            e_unp(:,j) = serr(sim.unps{j}{dimk}(t12,:),2);
            lr(:,j) = mean(sim.lrs{j}{dimk}(t12,:),2);            
            e_lr(:,j) = serr(sim.lrs{j}{dimk}(t12,:),2);
        end

        if dimk==2
            v = (sim.vals{j}{dimk}(N1,:));
            val1(:,j) = 1./(1+exp(-.5*v));
        end
        
        a{dimk}(:,j) = sim.lrs{j}{dimk}(ttend,:);
        
    end       
end

mv1 = mean(val1);


a = a{1};
mlr = mean(a,1);
elr = serr(a,1);

smax = 0.6;

labels = specs(2:3,:);
condition_names = labels(2,1:2);
xstr = {def('lr'), def('vol'), def('unp')};

fsy = def('fsy');

%----------------------
if nargin<1
    close all;    
    nr = 3;
    nc = 2;            
    fsiz = [0 0 .4 .7];    
    subplots = 1:6;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

h(1:2) = plot_bar(nr,nc,subplots(1:2),{m_emp',mlr'},{e_emp',elr'},labels,['Time in food cup %' xstr(1)]);
set(h(2),'ytick',0:.1:.6);
set(h(2),'ylim',[0 .6]);
title(h(1),'Data','fontsize',fsy);
title(h(2),'Model','fontsize',fsy);

[h(3:4),hp] = plot_signal(nr,nc,subplots(3:4),{vol(:,1:2),vol(:,3:4)},{e_vol(:,1:2),e_vol(:,3:4)},[xstr(2) xstr(2)],'',N1,[0 smax]);
legend(hp(1,:),condition_names,'fontsize',fsy,'location','northwest','box','off');

title(h(3),'Control');
title(h(4),'Lesioned');

h(5:6) = plot_signal(nr,nc,subplots(5:6),{unp(:,1:2),unp(:,3:4)},{e_unp(:,1:2),e_unp(:,3:4)},[xstr(3) xstr(3)],'',N1,[0 smax]);

%----------------------

tbl_data.empirical = ph1_emp;
tbl_data.simulated = mv1;
tbl_data.labels = labels;

q = [tbl_data.empirical' round(tbl_data.simulated'*10^4)/100];
end

%-----------------------
function [o, s]=timeseries(config,mode)
isshift = strcmpi(mode,'shift');

N1 = config.N1;
N2 = config.N2;
omega = config.omega;

x1 = zeros(N1,1);
j = randperm(N1);
j = j(1:(N1/2));
x1(j) = 1;

y = zeros(N2,1);
j = randperm(N2);
j = j(1:(N2/2));
y(j) = 1;

s1 = ones(size(x1)); 
x2 = ones(N2,1); 
s2 = ones(size(x2));
if ~isshift
    x2 = y;
else
    s2 = y;   
    x2 = ones(N2,1);
    x2(s2==0) = 0;
end

s = [s1; s2];
s = s + sqrt(omega)*randn(N1+N2,1);

x = [x1; x2];
o = x + sqrt(omega)*randn(N1+N2,1);
end
