function srf=constrSurfaceRevolve(crv,pnt_rot,vctr_rot,sang,eang)
% construct a surface by revolving a curve
% 
% input:
% crv(Curve): bus curve
% pnt_rot(vector): rotate center
% vctr_rot(vector): rotate vector
% theta(double): rotate angle(deg)
%
if nargin < 5
    eang=[];
    if nargin < 4
        sang=[];
        if nargin < 3
            vctr_rot=[];
            if nargin < 2
                pnt_rot=[];
            end
        end
    end
end

if isempty(eang), eang=2*pi;end
if isempty(sang), sang=0;end
if isempty(vctr_rot), vctr_rot=[1;0;0];end
if isempty(pnt_rot), pnt_rot=zeros(3,1);end

if length(pnt_rot) ~= 3
    error('constrSurfaceRevolve: all point and vector coordinates must be 3D');
end

if eang < sang
    temp=sang;
    sang=eang;
    eang=temp;
end

sweep=eang-sang; % sweep angle of arc
if sweep < 0
  sweep=2*pi+sweep;
end

vctr_rot=vctr_rot(:)./norm(vctr_rot);
if [0,0,1]*vctr_rot == 1 || [0,0,1]*vctr_rot == -1
    vctr_ref_x=[1,0,0];
else
    vctr_ref_x=cross([0;0;1],vctr_rot);
end
vctr_ref_y=cross(vctr_rot,vctr_ref_x);
rot_mat=[vctr_ref_x,vctr_ref_y,vctr_rot]'; % matrix to rotate curve into alignment with the z-axis

% translate and rotate the original curve into alignment with the z-axis
crv=crv.translate(-pnt_rot);
crv=crv.rotate(rot_mat);

% construct an arc 
arc=constrCurveCircular(1.0,[],0.0,sweep);
arc=arc.insertDimension(2,0);

% construct the revolved surface
coefs=zeros(arc.u_coef_num,crv.u_coef_num,4);
angle=atan2(crv.coefs(:,2),crv.coefs(:,1));angle(angle < 0)=2*pi+angle(angle < 0);
radius=vecnorm(crv.coefs(:,1:2),2,2);
for i=1:crv.u_coef_num
    coefs(:,i,:)=([radius(i),radius(i),0,1].*arc.coefs)*...
        tranH([0.0;0.0;crv.coefs(i,3)])'*rotZH(angle(i))';
    coefs(:,i,4)=coefs(:,i,4)*crv.coefs(i,4);
end
srf=Surface(coefs,{arc.u_knotvctr,crv.u_knotvctr});

% rotate and vectrans the surface back into position
srf=srf.rotate(rot_mat');
srf=srf.translate(pnt_rot);
srf=srf.rotate(axang2rotm([vctr_rot(:)',sang]),pnt_rot);

    function rot_mat=rotZH(angle)
        sn=sin(angle);
        cn=cos(angle);
        rot_mat=[cn -sn 0 0; sn cn 0 0; 0 0 1 0; 0 0 0 1];
    end

    function tran_mat=tranH(tran_vctr)
        tran_vctr=tran_vctr(:);
        tran_mat=[eye(length(tran_vctr)),tran_vctr;zeros(1,length(tran_vctr)),1];
    end
end
