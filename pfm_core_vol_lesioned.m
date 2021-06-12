function [val,vol,unp,lr] = pfm_core_vol_lesioned(o,x0_unc,lambda_u,v0,u0,nparticles)
dz = 0;
dy = length(lambda_u);

rng = bsxfun(@times,u0,ones(length(u0),2)).^-1;

state_model = @(particles)pf_state_transition(particles, lambda_u);
measurement_model = @pf_measurement_likelihood;

pf = particleFilter(state_model,measurement_model);
initialize(pf,nparticles,rng);

pf.StateEstimationMethod = 'mean';
pf.ResamplingMethod = 'systematic';

N = size(o,1);
D = size(o,2);

m = zeros(D,nparticles);
w = x0_unc*ones(D,nparticles);

estimated = nan(N,dz+dy);
val = nan(N,D);
lr = nan(N,D);
for t=1:N
    ot = o(t,:);
    
    estimated(t,:) = predict(pf);
    val(t,:) = pf.Weights*m';
    
    correct(pf,ot,v0,m,w);
    [m,w,k] = kalman(pf.Particles,ot,v0,m,w);
    lr(t,:) = pf.Weights*k';    
end

vol = bsxfun(@times,v0,ones(N,length(v0)));
unp = estimated(:,dz+(1:dy)).^-1;
end

%------------------------------
function particles = pf_state_transition(particles, lambdar)
% number of states x number of particles

% the mean of particles does not change for normal state-space model
lambda = lambdar;
eta = 1-lambda;
nu = .5./(1-eta);
for i=1:length(lambda)
    z = particles(i,:);
    epsil = betarnd(eta(i)*nu(i),(1-eta(i))*nu(i), size(z)) + eps;
    e = (eta(i).^-1)*epsil;    
    z = z.*e;
    particles(i,:) = z;
end
end

function likelihood = pf_measurement_likelihood(particles, measurement, v0, m, w)
% the mean of particles does not change for normal state-space model
D = length(measurement);

y = particles;
u = y.^-1;
v = bsxfun(@times,v0',ones(D,size(y,2)));
u = bsxfun(@times,u,ones(D,1));

l = nan(D,size(particles,2));
for i=1:length(measurement)
    l(i,:) = normpdf(measurement(i),m(i,:),sqrt(w(i,:)+v(i,:)+u(i,:)));
end

likelihood = prod(l,1);
end

function [m, w, k]=kalman(particles,o,v0,m,w)
D = length(o);
y = particles;
v = bsxfun(@times,v0',ones(D,size(y,2)));
u = y.^-1;

u = bsxfun(@times,u,ones(D,1));

k = (w+v)./(w+v+u);
m = m + k.*(o'-m);
w = (1-k).*(w+v);
end
