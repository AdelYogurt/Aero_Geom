function srf=constrSurface2(pnt_min,pnt_max,u_pole_num,v_pole_num)
% generate surface by diagonal line 2 tip point
%
if nargin < 4
    v_pole_num=[];
    if nargin < 3
        u_pole_num=[];
    end
end

bou=[pnt_min;pnt_max];
p11=[bou(1,1),bou(1,2)];
p21=[bou(2,1),bou(1,2)];
p12=[bou(1,1),bou(2,2)];
p22=[bou(2,1),bou(2,2)];

if isempty(u_pole_num), u_pole_num=2;end
if isempty(v_pole_num), v_pole_num=2;end

srf=constrSurface4(p11,p21,p12,p22,u_pole_num,v_pole_num);
end