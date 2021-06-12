function [val,vol,unp,lr] = pfm_core(o,x0_unc,lambda_v,lambda_u,v0,u0,nparticles)
dz = length(lambda_v);
dy = length(lambda_u);

vu0 = [v0 u0];
rng = bsxfun(@times,vu0',ones(length(vu0),2)).^-1;

state_model = @(particles)pf_state_transition(particles, lambda_v, lambda_u);
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
    correct(pf,ot,[dz, dy],m,w);
    [m,w,k] = kalman(pf.Particles,ot,[dz, dy],m,w);
    lr(t,:) = pf.Weights*k';
end

vol = estimated(:,1:dz).^-1;
unp = estimated(:,dz+(1:dy)).^-1;
end


%------------------------------
function particles = pf_state_transition(particles, lambdav, lambdar)
% number of states x number of particles

% the mean of particles does not change for normal state-space model
lambda = [lambdav lambdar];
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

function likelihood = pf_measurement_likelihood(particles, measurement, d, m, w)
% the mean of particles does not change for normal state-space model

z = particles(1:d(1),:);
y = particles(d(1)+(1:d(2)),:);
sigma = z.^-1;
omega = y.^-1;

D = length(measurement);

sigma = bsxfun(@times,sigma,ones(D,1));
omega = bsxfun(@times,omega,ones(D,1));

l = nan(D,size(z,2));
for i=1:length(measurement)
    l(i,:) = normpdf(measurement(i),m(i,:),sqrt(w(i,:)+sigma(i,:)+omega(i,:)));
end

likelihood = prod(l,1);
end

function [m, w, k]=kalman(particles,o,d,m,w)

z = particles(1:d(1),:);
y = particles(d(1)+(1:d(2)),:);
v = z.^-1;
u = y.^-1;

D = length(o);

v = bsxfun(@times,v,ones(D,1));
u = bsxfun(@times,u,ones(D,1));

k = (w+v)./(w+v+u);
m = m + k.*(o'-m);
w = (1-k).*(w+v);

end
