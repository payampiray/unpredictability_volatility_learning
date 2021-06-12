function [vol,unp,lr,val] = model_pfm(o,specs,islesioned)
if nargin<3, islesioned = 0; end

lambda_u = specs.lambda_u;
u = specs.u0;
np = specs.nparticles;

x0_unc = 1;
if isfield(specs,'x0_unc')
    x0_unc = specs.x0_unc;
end

if ~islesioned
    lambda_v = specs.lambda_v;
    v = specs.v0;
    [val,vol,unp,lr] = pfm_core(o,x0_unc,lambda_v,lambda_u,v,u,np);
else
    v = specs.v0_lesioned;
    [val,vol,unp,lr] = pfm_core_vol_lesioned(o,x0_unc,lambda_u,v,u,np);
end

dim = size(lr,2);
if size(vol,2)<dim
    vol = repmat(vol,1,dim);
end
if size(unp,2)<dim
    unp = repmat(unp,1,dim);
end

end
