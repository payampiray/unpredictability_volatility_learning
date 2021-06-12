function [val,vol,unp,lr,unc] = pf_core(o,x0_unc,lambda_v,lambda_u,v0,u0,nparticles)

z_rng = [v0 v0].^-1;
y_rng = [u0 u0].^-1;

state_model = @(particles)pf_state_transition(particles, lambda_v, lambda_u);
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

vol = estimated(:,1).^-1;
unp = estimated(:,2).^-1;
end


%------------------------------
function particles = pf_state_transition(particles, lambda_v, lambda_u)

z = particles(1,:);
eta = 1-lambda_v;
nu = .5/(1-eta);
epsil = betarnd(eta*nu,(1-eta)*nu, size(z)) + eps;
e = (eta.^-1)*epsil;
z = z.*e;

y = particles(2,:);
eta = 1-lambda_u;
nu = .5/(1-eta);
epsil = betarnd(eta*nu,(1-eta)*nu, size(y)) + eps;
e = (eta.^-1)*epsil;
y = y.*e;

particles = [z; y];

end

function likelihood = pf_measurement_likelihood(particles, measurement, m, w)

z = particles(1,:);
y = particles(2,:);
v = z.^-1;
u = y.^-1;

likelihood = normpdf(measurement,m,sqrt(w+v+u));

end

function [m, w, k]=kalman(particles,outcome,m,w)
z = particles(1,:);
y = particles(2,:);
v = z.^-1;
u = y.^-1;

k = (w+v)./(w+v+u);
m = m + k.*(outcome-m);
% w = (1-k).*(w+v);

w = u./(w+v+u).*(w+v);
end
