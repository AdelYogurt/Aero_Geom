function B_list=baseFcnB(us,k)
% Basis function for B-Spline
%
% input:
% u: parametric points
% k: spline degree
%
% output:
% B_list(matrix): Bezier-Basis functions vector, numel(u)*(k+1)
%
if any(us < 0) || any(us > 1)
    error('baseFcnB: u out of boundary of [0,1]')
end

if any(isnan(us))
    error('baseFcnB: u have nan')
end

u_num=numel(us);
B_list=zeros(u_num,k+1);
for i=0:k
    B_list(:,i+1)=nchoosek(k,i)*(us).^i.*(1-us).^(k-i);
end
end