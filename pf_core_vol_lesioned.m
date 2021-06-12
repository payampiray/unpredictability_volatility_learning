function [val,vol,unp,lr,unc] = pf_core_vol_lesioned(o,x0_unc,lambda_u,v0,u0,nparticles)
y_rng = [u0 u0].^-1;

state_model = @(particles)pf_state_transition(particles, lambda_u);
measurement_model = @pf_measurement_likelihood;

pf = particleFilter(state_model,measurement_model);
initialize(pf,nparticles,y_rng);

pf.StateEstimationMethod = 'mean';
pf.ResamplingMethod = 'systematic';

N = length(o);
estimated = nan(N,1);
val = nan(N,1);
unc = nan(N,1);
lr = nan(N,1);

m = zeros(1,nparticles);
w = x0_unc*ones(1,nparticles);

for t=1:size(o,1)    
    estimated(t) = predict(pf);
    correct(pf,o(t),m,w,v0);
    [m,w,k] = kalman(pf.Particles,o(t),m,w,v0);
    val(t) = pf.Weights*m';        
    unc(t) = pf.Weights*w';        
    lr(t) = pf.Weights*k';
end
unp = estimated(:,1).^-1;
vol = v0*ones(size(unp));
end


%------------------------------
function particles = pf_state_transition(particles, lambda_u)
% number of states x number of particles

% the mean of particles does not change for normal state-space model

y = particles(1,:);
eta = 1-lambda_u;
nu = .5/(1-eta);
epsil = betarnd(eta*nu,(1-eta)*nu, size(y)) + eps;
e = (eta.^-1)*epsil;
y = y.*e;

particles = y;

end

function likelihood = pf_measurement_likelihood(particles, measurement, m,w,v0)
% the mean of particles does not change for normal state-space model

y = particles(1,:);
v = v0;
u = y.^-1;

likelihood = normpdf(measurement,m,sqrt(w+v+u));

end

function [m, w, k]=kalman(particles,yt,m,w,v0)
y = particles(1,:);
v = v0;
u = y.^-1;

k = (w+v)./(w+v+u);
m = m + k.*(yt-m);
w = (1-k).*(w+v);
end
