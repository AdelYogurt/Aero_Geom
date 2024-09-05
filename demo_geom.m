clc;
clear;
close all hidden;

%% PATH

addpath src_constr\
addpath src_geom\

%% curve

group_num=5;deg=3;
pnts=[0,0,0;(0.5:0.5:group_num)',repmat([1;0],group_num,1),repmat([1;0],group_num,1)];
axe_hdl=axes(figure());
scatter3(axe_hdl,pnts(:,1),pnts(:,2),pnts(:,3),'Marker','*','MarkerEdgeColor',[0.9290 0.6940 0.1250]);

crv=Curve(pnts,deg);
crv.displayGeom(axe_hdl);
crv.displayPole(axe_hdl);

crv=interpPointToCurve(pnts,deg);
crv.displayGeom(axe_hdl,struct('color',[0.9290 0.6940 0.1250]));
axis equal;title('demo of curve');

%% CST curve

C_par=[0.5,1.0];S_par=[1,0.06];B_par=[0,0.0];
axe_hdl=axes(figure());

crv=CurveCST(C_par,S_par,B_par);
crv.displayGeom(axe_hdl,struct('LineStyle',':'));
axis equal;title('demo of CST curve');

pole_num=5;deg=4;
pnt_list=[linspace(0,1,pole_num)',rand(pole_num,1)*0.12];
crv=crv.addSpline(pnt_list,deg);
crv.displayGeom(axe_hdl);
crv.displayPole(axe_hdl);

crv=crv.convertSpline();
crv.displayGeom(axe_hdl,struct('Color','g','LineStyle','--'));
crv.displayPole(axe_hdl,struct('Color','k','LineStyle','-.','Marker','s'));

%% fit airfoil with CST curve

airfoil_data=importdata('src_geom/airfoil/NACA0012.txt');
airfoil_U=airfoil_data(1:66,:);airfoil_L=airfoil_data(67:end,:);

airfoil_data=importdata('src_geom/airfoil/NACA4412.txt');
airfoil_U=airfoil_data(1:18,:);airfoil_L=airfoil_data(19:end,:);

airfoil_data=importdata('src_geom/airfoil/RAE2822.txt');
airfoil_U=airfoil_data(1:65,:);airfoil_L=airfoil_data(66:end,:);

airfoil_data=importdata('src_geom/airfoil/Clark_Y.txt');
airfoil_U=airfoil_data(1:61,:);airfoil_L=airfoil_data(62:end,:);

LX=1;LY=0.06;C_par=[0.5,1];
pole_num=5;deg=3;

crv_L=CurveCST(C_par,[LX,-LY],[0,airfoil_L(end,2)]);
crv_U=CurveCST(C_par,[LX,LY],[0,airfoil_U(end,2)]);

crv_L=crv_L.fitSpline(airfoil_L,deg,pole_num,airfoil_L(:,1));
crv_U=crv_U.fitSpline(airfoil_U,deg,pole_num,airfoil_U(:,1));
disp(['fit_err_L: ',num2str(crv_L.fitError)]);
disp(['fit_err_U: ',num2str(crv_U.fitError)]);

[fit_err_L,crv_L]=crv_L.optimClass();
[fit_err_U,crv_U]=crv_U.optimClass();
disp(['fit_err_L: ',num2str(fit_err_L)]);
disp(['fit_err_U: ',num2str(fit_err_U)]);

axe_hdl=axes(figure());
line(axe_hdl,airfoil_L(:,1),airfoil_L(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.8500 0.3250 0.0980])
line(axe_hdl,airfoil_U(:,1),airfoil_U(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.8500 0.3250 0.0980])
crv_L.displayGeom(axe_hdl);
crv_U.displayGeom(axe_hdl);
axis equal;title('demo of CST curve fitting airfoil');

writematrix(crv_L.getPoles(),'CSTshape_L.txt');
writematrix(crv_U.getPoles(),'CSTshape_U.txt');
writematrix(crv_L.u_knotvctr,'CSTshape_L_knotvctr.txt');
writematrix(crv_U.u_knotvctr,'CSTshape_U_knotvctr.txt');

%% surface

pnt_num=5;u_deg=3;v_deg=3;
[pnts_x,pnts_y]=ndgrid(linspace(0,1,pnt_num),linspace(0,1,pnt_num));pnts_z=round(pnts_x+pnts_y)*0.1;
axe_hdl=axes(figure());
surface(axe_hdl,pnts_x,pnts_y,pnts_z,'Marker','*','MarkerEdgeColor','r','LineStyle','none','FaceAlpha',0);
pnts=cat(3,pnts_x,pnts_y,pnts_z);

srf=Surface(pnts,u_deg,v_deg);
srf.displayGeom(axe_hdl);
srf.displayPole(axe_hdl);
axis equal;view(3);title('demo of surface');

srf=interpPointToSurface(pnts,u_deg,v_deg);
srf.displayGeom(axe_hdl);

%% CST surface

C_par_x=[0,0];C_par_y=[0.5,0];C_par_zv=[15,15];C_par_zu=[0.5,0];
S_par=[2,0.5,0.3];B_par=[0,0.0];sym_x=false;sym_y=true;
axe_hdl=axes(figure());
srf=SurfaceCST(C_par_x,C_par_y,C_par_zv,C_par_zu,S_par,B_par,sym_x,sym_y);
srf.displayGeom(axe_hdl,struct('LineStyle','none','FaceAlpha',0.2));
axis equal;view(3);title('demo of CST surface');

pnt_num=5;u_deg=3;v_deg=3;
[pnts_x,pnts_y]=ndgrid(linspace(0,1,pnt_num)*2,linspace(0,1,pnt_num)*0.5);pnts_z=rands(pnt_num,pnt_num)*0.2+0.2;
pnts=cat(3,pnts_x,pnts_y,pnts_z);
srf=srf.addSpline(pnts,u_deg,v_deg);
srf.displayGeom(axe_hdl);
srf.displayPole(axe_hdl);

srf=srf.convertSpline();
srf.displayGeom(axe_hdl,struct('FaceColor','g','FaceAlpha',0.5,'LineStyle','none'));
srf.displayPole(axe_hdl,struct('FaceAlpha',0,'MarkerEdgeColor','k','LineStyle','-.','Marker','s'));

%% mapping generate surface

axe_hdl=axes(figure());

line_u1=[0,0,0;0,1,0;0.5,2,0];
line_u2=[0,0,2.5;0.5,1,2;1,2.5,2];
line_1v=[0,0,0;0,1,1;0,0,2.5];
line_2v=[0.5,2,0;0.5,1.5,1;1,2.5,2];
crv_u1=Curve(line_u1);
crv_u2=Curve(line_u2);
crv_1v=Curve(line_1v);
crv_2v=Curve(line_2v);
crv_u1.displayGeom(axe_hdl);
crv_u2.displayGeom(axe_hdl);
crv_1v.displayGeom(axe_hdl);
crv_2v.displayGeom(axe_hdl);

srf=constrSurfaceBound(crv_u1,crv_u2,crv_1v,crv_2v);

srf.displayGeom(axe_hdl,struct('LineStyle','none','FaceAlpha',0.5));
srf.displayPole(axe_hdl);
axis equal;view(3);title('demo of bounding form surface');

%% fit wing with CST surface

load('src_geom/wing/wing_RAE2822_NACA4412.mat');
u_degree=3;v_degree=3;u_pole_num=4;v_pole_num=8;
C_par_x=[0,0];C_par_y=[0,0];C_par_zv=[0,0];C_par_zu=[0.5,1];
srf_U=SurfaceCST(C_par_x,C_par_y,C_par_zv,C_par_zu);
srf_L=SurfaceCST(C_par_x,C_par_y,C_par_zv,C_par_zu);

u_nodes_U=X_U(:,1)/max(X_U(:,1));
v_nodes_U=Y_U(1,:)/max(Y_U(1,:));
u_nodes_L=X_L(:,1)/max(X_L(:,1));
v_nodes_L=Y_L(1,:)/max(Y_L(1,:));

srf_U=srf_U.fitSpline(cat(3,X_U,Y_U,Z_U),u_degree,v_degree,u_pole_num,v_pole_num,u_nodes_U,v_nodes_U);
srf_L=srf_L.fitSpline(cat(3,X_L,Y_L,Z_L),u_degree,v_degree,u_pole_num,v_pole_num,u_nodes_L,v_nodes_L);

axe_hdl=axes(figure());
srf_U.displayGeom(axe_hdl);
srf_U.displayPole(axe_hdl);
srf_L.displayGeom(axe_hdl);
srf_L.displayPole(axe_hdl);
axis equal;view(3);title('demo of CST surface fitting wing');

%% volume

pnt_num=5;u_deg=3;v_deg=3;w_deg=1;
[pnts_x,pnts_y]=ndgrid(linspace(0,1,pnt_num),linspace(0,1,pnt_num));
pnts_z1=round(pnts_x+pnts_y)*0.1;pnts_z2=-round(pnts_x+pnts_y)*0.1+1.0;
axe_hdl=axes(figure());
surface(axe_hdl,pnts_x,pnts_y,pnts_z1,'Marker','*','MarkerEdgeColor','r','LineStyle','none','FaceAlpha',0);
surface(axe_hdl,pnts_x,pnts_y,pnts_z2,'Marker','*','MarkerEdgeColor','r','LineStyle','none','FaceAlpha',0);
pnts=cat(4,repmat(pnts_x,[1,1,2]),repmat(pnts_y,[1,1,2]),cat(3,pnts_z1,pnts_z2));

vol=Volume(pnts,u_deg,v_deg,w_deg);
vol.displayGeom(axe_hdl);
vol.displayPole(axe_hdl);
view(3);axis equal;title('demo of volume');
