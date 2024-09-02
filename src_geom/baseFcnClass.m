function [c,dc_du]=baseFcnClass(u,N1,N2)
% default of class function
%
% input:
% u: parametric points
% N1: class parameter 1
% N2: class parameter 2
%
% output:
% c: class value
% dc_du: gradient value
%
NC=calNormPar(N1,N2);
u_N1=u.^N1;u_N2=(1-u).^N2;
c=u_N1.*u_N2./NC;

if nargout > 1
    dc_du=(N1.*u.^(N1-1).*u_N2-N2.*u_N1.*(1-u).^(N2-1))./NC;
end
end

function nomlz_par=calNormPar(N1,N2)
% calculate normailize class function parameter by N1, N2
%
nomlz_par=(N1./(N1+N2)).^N1.*(N2./(N1+N2)).^N2;
nomlz_par((N1 == 0) & (N2 == 0))=1;
end
