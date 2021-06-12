function [val,vol,unp,lr,unc] = pf_core_norm(o,x0_unc,sigma_v,sigma_u,v0,u0,nparticles)


z_rng = log([v0 v0]);
y_rng = log([u0 u0]);

state_model = @(particles)pf_state_transition(particles, sigma_v, sigma_u);
measurement_model = @pf_measurement_likelihood;

pf = particleFilter(state_model,measurement_model);
initialize(pf,nparticles,[z_rng; y_rng]);

pf.StateEstimationMethod = 'mean';
pf.ResamplingMethod = 'systematic';

N = length(o);
estimated = nan(N,2);
val = nan(N,1);
unc = nan(N,1);
lr = nan(N,1);


m = zeros(1,nparticles);
w = x0_unc*ones(1,nparticles);

for t=1:size(o,1)    
    estimated(t,:) = predict(pf);    
    correct(pf,o(t),m,w);
    [m,w,k] = kalman(pf.Particles,o(t),m,w);
    val(t) = pf.Weights*m';
    unc(t) = pf.Weights*w';
    lr(t) = pf.Weights*k';
end

vol = exp(estimated(:,1));
unp = exp(estimated(:,2));
end


%------------------------------
function particles = pf_state_transition(particles, sigma_v, sigma_u)

v = particles(1,:);
epsil = randn(size(v));
v = v + sigma_v*epsil;

u = particles(2,:);
epsil = randn(size(u));
u = u + sigma_u*epsil;

particles = [v; u];

end

function likelihood = pf_measurement_likelihood(particles, measurement, m, w)

v = exp(particles(1,:));
u = exp(particles(2,:));

likelihood = normpdf(measurement,m,sqrt(w+v+u));

end

function [m, w, k]=kalman(particles,outcome,m,w)
v = exp(particles(1,:));
u = exp(particles(2,:));

k = (w+v)./(w+v+u);
m = m + k.*(outcome-m);
w = (1-k).*(w+v);
end
