function vol=constrVolume2(pnt_min,pnt_max,u_pole_num,v_pole_num,w_pole_num)
% generate volume by diagonal line 2 tip point
%
if nargin < 5
    w_pole_num=[];
    if nargin < 4
        v_pole_num=[];
        if nargin < 3
            u_pole_num=[];
        end
    end
end

bou=[pnt_min;pnt_max];
p111=[bou(1,1),bou(1,2),bou(1,3)];
p211=[bou(2,1),bou(1,2),bou(1,3)];
p121=[bou(1,1),bou(2,2),bou(1,3)];
p221=[bou(2,1),bou(2,2),bou(1,3)];
p112=[bou(1,1),bou(1,2),bou(2,3)];
p212=[bou(2,1),bou(1,2),bou(2,3)];
p122=[bou(1,1),bou(2,2),bou(2,3)];
p222=[bou(2,1),bou(2,2),bou(2,3)];

if isempty(u_pole_num), u_pole_num=2;end
if isempty(v_pole_num), v_pole_num=2;end
if isempty(w_pole_num), w_pole_num=2;end

vol=constrVolume8(p111,p211,p121,p221,p112,p212,p122,p222,u_pole_num,v_pole_num,w_pole_num);
end