function srf=constrSurfaceLoft(crv_u1,crv_u2,crv_1v,crv_2v,geom_torl)
% curve to generate surface by solving poisson equation
%
if nargin < 5
    geom_torl=[];
    if nargin < 4
        crv_2v=[];
        if nargin < 3
            crv_1v=[];
        end
    end
end
if isempty(geom_torl), geom_torl=1e-6;end

if (~isempty(crv_u1) && ~isempty(crv_u2)) || (~isempty(crv_1v) && ~isempty(crv_2v))
    if ~isempty(crv_u1) && ~isempty(crv_u2)
        [crv_u1,crv_u2,u_order,u_knotvctr]=curveMatch(crv_u1,crv_u2);
        coefs_u1=crv_u1.coefs;coefs_u2=crv_u2.coefs;
    end

    if ~isempty(crv_1v) && ~isempty(crv_2v)
        [crv_1v,crv_2v,v_order,v_knotvctr]=curveMatch(crv_1v,crv_2v);
        coefs_1v=crv_1v.coefs;coefs_2v=crv_2v.coefs;
    end

    if isempty(crv_1v) && isempty(crv_2v)
        coefs_1v=[crv_u1.coefs(1,:);crv_u2.coefs(1,:)];
        coefs_2v=[crv_u1.coefs(end,:);crv_u2.coefs(end,:)];
        v_order=1;v_knotvctr=[0,0,1,1];
    elseif isempty(crv_u1) && isempty(crv_u2)
        coefs_u1=[crv_1v.coefs(1,:);crv_2v.coefs(1,:)];
        coefs_u2=[crv_1v.coefs(end,:);crv_2v.coefs(end,:)];
        u_order=1;u_knotvctr=[0,0,1,1];
    elseif isempty(crv_1v) || isempty(crv_2v)
        crv_v=[crv_1v,crv_2v];
        v_order=crv_v.u_order;
        v_knotvctr=crv_v.u_knotvctr;

        if ~isempty(crv_1v)
            % generate new ctrl
            crv_2v=Curve([coefs_u1(end,:);coefs_u2(end,:)],[0,0,1,1]);
            [~,crv_2v]=curveMatch(crv_v,crv_2v);
            coefs_1v=crv_v.coefs;
            coefs_2v=crv_2v.coefs;
        elseif ~isempty(crv_2v)
            % generate new ctrl
            crv_1v=Curve([coefs_u1(1,:);coefs_u2(1,:)],[0,0,1,1]);
            [~,crv_1v]=curveMatch(crv_v,crv_1v);
            coefs_1v=crv_1v.coefs;
            coefs_2v=crv_v.coefs;
        end
    elseif isempty(crv_u1) || isempty(crv_u2)
        crv_u=[crv_u1,crv_u2];
        u_order=crv_u.u_order;
        u_knotvctr=crv_u.u_knotvctr;

        if ~isempty(crv_u1)
            % generate new ctrl
            crv_u2=Curve([coefs_1v(end,:);coefs_1v(end,:)],[0,0,1,1]);
            [~,crv_u2]=curveMatch(crv_u,crv_u2);
            coefs_u1=crv_u.coefs;
            coefs_u2=crv_u2.coefs;
        elseif ~isempty(crv_u2)
            % generate new ctrl
            crv_u1=Curve([coefs_1v(1,:);coefs_1v(1,:)],[0,0,1,1]);
            [~,crv_u1]=curveMatch(crv_u,crv_u1);
            coefs_u1=crv_u1.coefs;
            coefs_u2=crv_u.coefs;
        end
    end
else
    error('constrSurfaceBound: error input format');
end

u_pole_num=size(coefs_u1,1);v_pole_num=size(coefs_1v,1);
u_list=interp1(linspace(0,1,u_pole_num-u_order+1),u_knotvctr(u_order+1:u_pole_num+1),linspace(0,1,u_pole_num));
v_list=interp1(linspace(0,1,v_pole_num-v_order+1),v_knotvctr(v_order+1:v_pole_num+1),linspace(0,1,v_pole_num));
[U2,V2]=ndgrid(u_list,v_list);
U1=1-U2;V1=1-V2;
% Ctrl=MapGrid(Ctrl_u1,Ctrl_u2,Ctrl_1v,Ctrl_2v,u_list,v_list);

% generata U, V and UV direction coefs
coefs_u1=reshape(coefs_u1,size(coefs_u1,1),1,size(coefs_u1,2));
coefs_u2=reshape(coefs_u2,size(coefs_u2,1),1,size(coefs_u2,2));
coefs_1v=reshape(coefs_1v,1,size(coefs_1v,1),size(coefs_1v,2));
coefs_2v=reshape(coefs_2v,1,size(coefs_2v,1),size(coefs_2v,2));
ctrl_11=coefs_u1(1,1,:);ctrl_21=coefs_u1(end,1,:);
ctrl_12=coefs_u2(1,1,:);ctrl_22=coefs_u2(end,1,:);
coefs_U=coefs_u1.*V1+coefs_u2.*V2;
coefs_V=coefs_1v.*U1+coefs_2v.*U2;
coefs_UV=...
    ctrl_11.*V1.*U1+ctrl_21.*V1.*U2+...
    ctrl_12.*V2.*U1+ctrl_22.*V2.*U2;

% combine coefficient to construct Coons surface
coefs=coefs_U+coefs_V-coefs_UV;

srf=Surface(coefs,{u_knotvctr,v_knotvctr});
end
