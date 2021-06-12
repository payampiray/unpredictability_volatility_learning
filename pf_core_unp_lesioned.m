function [val,vol,unp,lr,unc] = pf_core_unp_lesioned(o,x0_unc,lambda_v,v0,u0,nparticles)
z_rng = [v0 v0].^-1;

state_model = @(particles)pf_state_transition(particles, lambda_v);
measurement_model = @pf_measurement_likelihood;

pf = particleFilter(state_model,measurement_model);
initialize(pf,nparticles,z_rng);

pf.StateEstimationMethod = 'mean';
pf.ResamplingMethod = 'systematic';

N = length(o);
estimated = nan(N,1);
val = nan(N,1);
lr = nan(N,1);
unc = nan(N,1);

m = zeros(1,nparticles);
w = x0_unc*ones(1,nparticles);

for t=1:size(o,1)    
    estimated(t) = predict(pf);
    correct(pf,o(t),m,w,u0);
    [m,w,k] = kalman(pf.Particles,o(t),m,w,u0);
    val(t) = pf.Weights*m';
    unc(t) = pf.Weights*w';
    lr(t) = pf.Weights*k';
end
vol = estimated(:,1).^-1;
unp = u0*ones(size(vol));
end


%------------------------------
function particles = pf_state_transition(particles, lambda_v)
% number of states x number of particles

% the mean of particles does not change for normal state-space model

z = particles(1,:);
eta = 1-lambda_v;
nu = .5/(1-eta);
epsil = betarnd(eta*nu,(1-eta)*nu, size(z)) + eps;
e = (eta.^-1)*epsil;
z = z.*e;

particles = z;

end

function likelihood = pf_measurement_likelihood(particles, measurement, m,w,u0)
% the mean of particles does not change for normal state-space model

z = particles(1,:);
u = u0;
v = z.^-1;

likelihood = normpdf(measurement,m,sqrt(w+v+u));

end

function [m, w, k]=kalman(particles,yt,m,w,u0)
z = particles(1,:);
u = u0;
v = z.^-1;

k = (w+v)./(w+v+u);
m = m + k.*(yt-m);
w = (1-k).*(w+v);
end
