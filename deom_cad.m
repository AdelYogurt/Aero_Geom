clc;
clear;
close all hidden;

%% PATH

addpath src_cad\
addpath src_geom\
addpath src_topo\

%% read geom from IGES

geom_list=readGeomIGES('DPW4.iges');

for geom_idx=1:length(geom_list)
    geom=geom_list{geom_idx};
    geom.displayGeom();
end
view(3);axis equal;

%% write geom from IGES

atllas=ATLLAS();
atllas.displayGeom([],[],[],2);
atllas.displayDirect();
view(3);axis equal;

writeGeomIGES(num2cell(atllas.srf_list),'atllas.iges');

%% read geom from STEP

geom_list=readGeomIGES('DPW4.step');

for geom_idx=1:length(geom_list)
    geom=geom_list{geom_idx};
    geom.displayGeom();
end
view(3);axis equal;

%% write geom from STEP

atllas=ATLLAS();
atllas.displayGeom([],[],[],2);
atllas.displayDirect();
view(3);axis equal;

writeGeomIGES(num2cell(atllas.srf_list),'atllas.step');
