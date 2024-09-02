function srf=constrSurface4(p11,p21,p12,p22,u_pole_num,v_pole_num)
% generate surface by 4 tip point
%
if nargin < 6
    v_pole_num=[];
    if nargin < 5
        u_pole_num=[];
    end
end

p11=reshape(p11,1,1,[]);
p21=reshape(p21,1,1,[]);
p12=reshape(p12,1,1,[]);
p22=reshape(p22,1,1,[]);

if isempty(u_pole_num), u_pole_num=2;end
if isempty(v_pole_num), v_pole_num=2;end

[U2,V2]=ndgrid(linspace(0,1,u_pole_num),linspace(0,1,v_pole_num));
U2=reshape(U2,[size(U2),1]);
V2=reshape(V2,[size(V2),1]);
U1=1-U2;V1=1-V2;

pnts=...
    U1.*V1.*p11+...
    U2.*V1.*p21+...
    U1.*V2.*p12+...
    U2.*V2.*p22;

srf=Surface(pnts);
end