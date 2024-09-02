function crv=constrCurve2(p1,p2,u_pole_num)
% generate surface by 2 tip point
%
if nargin < 3
    u_pole_num=[];
end

p1=reshape(p1,1,[]);
p2=reshape(p2,1,[]);

if isempty(u_pole_num), u_pole_num=2;end

U2=linspace(0,1,u_pole_num);
U1=1-U2;

pnts=...
    U1.*p1+...
    U2.*p2;

crv=Curve(pnts);
end