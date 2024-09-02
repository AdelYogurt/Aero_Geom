function crv=constrCurveCircular(radius,cntr,sang,eang)
% construct a circular arc.
% 
% input:
% radius(optional): radius of the circle, default 1.0
% cntr(optional): center of the circle, default (0,0,0)
% sang(optional): start angle, default 0 radians (0 rad)
% eang(optional): end angle, default 2*pi radians (2*pi rad)
%
if nargin < 1, radius=1;end
if nargin < 2, cntr=[];end
if nargin < 4
  sang=0;
  eang=2*pi;
end

if eang < sang
    temp=sang;
    sang=eang;
    eang=temp;
end

sweep=eang-sang; % sweep angle of arc
sangd=sang;
if sweep < 0
  sweep=2*pi+sweep;
end
     
if abs(sweep) <= pi/2
  narcs=1; % number of arc segments
  u_knotvctr=[0 0 0 1 1 1];
elseif abs(sweep) <= pi
  narcs=2;
  u_knotvctr=[0 0 0 0.5 0.5 1 1 1];
elseif abs(sweep) <= 3*pi/2
  narcs=3;
  u_knotvctr=[0 0 0 1/3 1/3 2/3 2/3 1 1 1];
else
  narcs=4;
  u_knotvctr=[0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1];
end

dsweep=sweep/(2*narcs);     % arc segment sweep angle/2

% determine middle control point and weight
wm=cos(dsweep);
x =radius*wm;
y =radius*sin(dsweep);
xm=x+y*tan(dsweep);

% arc segment control points
% [w*x,w*y,w]
ctrlpt=[
    x,    -y,  1;
    wm*xm, 0, wm;
    x,     y,  1;];
     
% build up complete arc from rotated segments
coefs=zeros(2*narcs+1,3);   % nurbs control points of arc
xx=rotZ(sangd+dsweep);
coefs(1:3,:)=ctrlpt*xx';     % rotate to start angle
xx=rotZ(2*dsweep);
for n=2:narcs
   m=2*n+[0 1];
   coefs(m,:)=coefs(m-2,:)*xx';
end

crv=Curve(coefs,u_knotvctr);

% vectrans arc if necessary
if ~isempty(cntr) 
  crv=crv.translate(cntr);
end

    function rotm=rotZ(angle)
        sn=sin(angle);
        cn=cos(angle);
        rotm=[cn -sn 0; sn cn 0; 0 0 1];
    end
end

% % demo
% for r=1:9
%   crv=constrCurveCircular(r,[r,0,0],45*pi/180,315*pi/180);
%   crv.displayGeom()
%   hold on;
% end
% hold off;
% axis equal;
% title('NURBS construction of several 2D arcs.');
