function [vol,unp,lr,val,unc] = model_pf(o,specs,islesioned)
if nargin<3, islesioned = 'healthy'; end

np = specs.nparticles;

x0_unc = 1;
if isfield(specs,'x0_unc')
    x0_unc = specs.x0_unc;
end    

switch lower(islesioned(1:3))
    case {'hea'}
        lambda_v = specs.lambda_v;
        lambda_u = specs.lambda_u;        
        v = specs.v0;
        u = specs.u0;
        [val,vol,unp,lr,unc] = pf_core(o,x0_unc,lambda_v,lambda_u,v,u,np);        
    case 'vol'
        lambda_u = specs.lambda_u;
        v_lesioned = specs.v0_lesioned;        
        u = specs.u0;
        [val,vol,unp,lr,unc] = pf_core_vol_lesioned(o,x0_unc,lambda_u,v_lesioned,u,np);        
    case 'unp'
        lambda_v = specs.lambda_v;
        u_lesioned = specs.u0_lesioned;
        v = specs.v0;
        [val,vol,unp,lr,unc] = pf_core_unp_lesioned(o,x0_unc,lambda_v,v,u_lesioned,np);        
    otherwise
        error('bad 3rd input: %d',islesioned);
end


end