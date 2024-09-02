function vol=constrVolume8(p111,p211,p121,p221,p112,p212,p122,p222,...
    u_pole_num,v_pole_num,w_pole_num)
% generate volume by 8 tip point
%
if nargin < 11
    w_pole_num=[];
    if nargin < 10
        v_pole_num=[];
        if nargin < 9
            u_pole_num=[];
        end
    end
end

p111=reshape(p111,1,1,1,[]);
p211=reshape(p211,1,1,1,[]);
p121=reshape(p121,1,1,1,[]);
p221=reshape(p221,1,1,1,[]);
p112=reshape(p112,1,1,1,[]);
p212=reshape(p212,1,1,1,[]);
p122=reshape(p122,1,1,1,[]);
p222=reshape(p222,1,1,1,[]);

if isempty(u_pole_num), u_pole_num=2;end
if isempty(v_pole_num), v_pole_num=2;end
if isempty(w_pole_num), w_pole_num=2;end

[U2,V2,W2]=ndgrid(linspace(0,1,u_pole_num),linspace(0,1,v_pole_num),linspace(0,1,w_pole_num));
U2=reshape(U2,[size(U2),1]);
V2=reshape(V2,[size(V2),1]);
W2=reshape(W2,[size(W2),1]);
U1=1-U2;V1=1-V2;W1=1-W2;

pnts=...
    U1.*V1.*W1.*p111+...
    U2.*V1.*W1.*p211+...
    U1.*V2.*W1.*p121+...
    U2.*V2.*W1.*p221+...
    U1.*V1.*W2.*p112+...
    U2.*V1.*W2.*p212+...
    U1.*V2.*W2.*p122+...
    U2.*V2.*W2.*p222;

vol=Volume(pnts);
end