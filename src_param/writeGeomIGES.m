function writeGeomIGES(geom_list,iges_filestr)
% write geometry entity into IGES file
%
if nargin < 2,iges_filestr=[];end
if isempty(iges_filestr),iges_filestr='geom.iges';end

if ~iscell(geom_list),geom_list={geom_list};end

[~,iges_filename,~]=fileparts(iges_filestr);
geom_num=length(geom_list);

iges_file=fopen(iges_filestr,'w');

% write the head statement
fprintf(iges_file,"                                                                        S      1\n");


fprintf(iges_file,"1H,,1H;,11H    1,38HD,32,308,15,308,15,11H    1,                        G      1\n");
fprintf(iges_file,"1.,2,2HMM,50,0.125,1E-08,499990.,,11,0;                                 G      2\n");

D_num=1;
P_num=1;

for geom_idx=1:geom_num
    geom=geom_list{geom_idx};

    if isa(geom,'Curve')
        [P_num,D_num]=writeCurveDir(geom,iges_file,D_num,P_num);
    elseif isa(geom,'CurveCST')
        geom=geom.convertSpline();
        [P_num,D_num]=writeCurveDir(geom,iges_file,D_num,P_num);
    elseif isa(geom,'Surface')
        [P_num,D_num]=writeSurfaceDir(geom,iges_file,D_num,P_num);
    elseif isa(geom,'SurfaceCST')
        geom=geom.convertSpline();
        [P_num,D_num]=writeSurfaceDir(geom,iges_file,D_num,P_num);
    elseif isa(geom,'Shape')
        % shape is a container of curve and surface
        crv_list=geom.crv_list;
        srf_list=geom.srf_list;
    else
        error('writeGeomIGES: unsupport geometry format');
    end
end

P_num=1;
counter=1;

for geom_idx=1:geom_num
    geom=geom_list{geom_idx};

    if isa(geom,'Curve')
        [P_num,counter]=writeCurvePar(geom,iges_file,P_num,counter);
    elseif isa(geom,'CurveCST')
        geom=geom.convertSpline();
        [P_num,counter]=writeCurvePar(geom,iges_file,P_num,counter);
    elseif isa(geom,'Surface')
        [P_num,counter]=writeSurfacePar(geom,iges_file,P_num,counter);
    elseif isa(geom,'SurfaceCST')
        geom=geom.convertSpline();
        [P_num,counter]=writeSurfacePar(geom,iges_file,P_num,counter);
    elseif isa(geom,'Shape')
        % shape is a container of curve and surface
        crv_list=geom.crv_list;
        srf_list=geom.srf_list;
    else
        error('writeGeomIGES: unsupport geometry format');
    end
end

% write the end statement
fprintf(iges_file,"S%7dG%7dD%7dP%7d%40sT%6s1\n",1,4,D_num-1,counter-1," "," ");
fclose(iges_file);
clear('iges_file');
end

%% geometry IGES

function [P_num,D_num]=writeCurveDir(crv,iges_file,D_num,P_num,FLAG_2D)
% write curve head information into IGES file
%
if nargin < 5,FLAG_2D=[];end
if isempty(FLAG_2D),FLAG_2D=false();end

if crv.coef_dim ~= 3 && crv.coef_dim ~= 4
    error('writeCurveDir: must have 2 or 3 dimensions to write to IGES file');
end

param_geom=6+length(crv.u_knotvctr)+crv.u_coef_num+3*crv.u_coef_num+5;

param_num=floor((param_geom-11)/3)+2;
if mod(param_geom-11,3) ~= 0
    param_num=param_num+1;
end

if FLAG_2D
    fprintf(iges_file,'     126%8d       0       0       1       0       0       001010501D%7d\n',P_num,D_num);
    fprintf(iges_file,'     126       0       0%8d       0                               0D%7d\n',param_num,D_num+1);
else
    fprintf(iges_file,'     126%8d       0       0       1       0       0       000000001D%7d\n',P_num,D_num);
    fprintf(iges_file,'     126       0       0%8d       0                               0D%7d\n',param_num,D_num+1);
end

D_num=D_num+2;
P_num=P_num+param_num;
end

function [P_num,counter]=writeCurvePar(crv,iges_file,P_num,counter)
% write curve parameter information info IGES file
%

[poles,weights]=crv.getPoles();
if all(weights == 1.0),bool_poly=1;
else,bool_poly=0;end

fprintf(iges_file,'%10d,%10d,%10d,%5d,%5d,%5d,%5d,        %7dP%7d\n',...
    126,crv.u_coef_num-1,crv.u_order,0,0,bool_poly,0,P_num,counter);
counter=counter+1;
pos_counter=0;

for knots_idx=1:length(crv.u_knotvctr)
    pos_counter=pos_counter+1;
    fprintf(iges_file,'%20.12g,',real(crv.u_knotvctr(knots_idx)));
    if mod(pos_counter,3) == 0
        fprintf(iges_file,'  %7dP%7d\n',P_num,counter);
        counter=counter+1;
        pos_counter=0;
    end
end

for u_idx=1:crv.u_coef_num
    pos_counter=pos_counter+1;
    fprintf(iges_file,'%20.12g,',real(weights(u_idx,end)));
    if mod(pos_counter,3) == 0
        fprintf(iges_file,'  %7dP%7d\n',P_num,counter);
        counter=counter+1;
        pos_counter=0;
    end
end

for u_idx=1:crv.u_coef_num
    for dim_idx=1:3
        pos_counter=pos_counter+1;
        fprintf(iges_file,'%20.12g,',real(poles(u_idx,dim_idx)));
        if mod(pos_counter,3) == 0
            fprintf(iges_file,'  %7dP%7d\n',P_num,counter);
            counter=counter+1;
            pos_counter=0;
        end
    end
end

if pos_counter == 1
    fprintf(iges_file,'%s    %7dP%7d\n',repmat(' ',1,40),P_num,counter);
    counter=counter+1;
elseif pos_counter == 2
    fprintf(iges_file,'%s    %7dP%7d\n',repmat(' ',1,20),P_num,counter);
    counter=counter+1;
end

fprintf(iges_file,'%20.12g,%20.12g,0.0,0.0,0.0;         ',min(crv.u_knotvctr),max(crv.u_knotvctr));
fprintf(iges_file,'  %7dP%7d\n',P_num,counter);
counter=counter+1;
P_num=P_num+2;
end

function [P_num,D_num]=writeSurfaceDir(srf,iges_file,D_num,P_num)
% write surface head information into IGES file
%

% A simpler calc based on cmlib definitions The 13 is for the
% 9 parameters at the start, and 4 at the end. See the IGES
% 5.3 Manual paraEntries=13+Knotsu+Knotsv+Weights+control points
if srf.coef_dim ~= 4
    error('Must have 3 dimensions to write to IGES file');
end
param_geom=13+length(srf.u_knotvctr)+length(srf.v_knotvctr)+srf.u_coef_num*srf.v_coef_num+3*srf.u_coef_num*srf.v_coef_num+1;

param_num=floor((param_geom-10)/3)+2;
if mod(param_geom-10,3) ~= 0
    param_num=param_num+1;
end

fprintf(iges_file,'     128%8d       0       0       1       0       0       000000001D%7d\n',P_num,D_num);
fprintf(iges_file,'     128       0       0%8d       0                               0D%7d\n',param_num,D_num+1);
D_num=D_num+2;
P_num=P_num+param_num;
end

function [P_num,counter]=writeSurfacePar(srf,iges_file,P_num,counter)
% write surface parameter information into IGES file
%
fprintf(iges_file,'%10d,%10d,%10d,%10d,%10d,          %7dP%7d\n',128,srf.u_coef_num-1,srf.v_coef_num-1,srf.u_order,srf.v_order,P_num,counter);
counter=counter+1;

[poles,weights]=srf.getPoles();
if all(weights == 1.0),bool_poly=1;
else,bool_poly=0;end

fprintf(iges_file,'%10d,%10d,%10d,%10d,%10d,          %7dP%7d\n',0,0,bool_poly,0,0,P_num,counter);

counter=counter+1;
pos_counter=0;

for u_idx=1:length(srf.u_knotvctr)
    pos_counter=pos_counter+1;
    fprintf(iges_file,'%20.12g,',real(srf.u_knotvctr(u_idx)));
    if mod(pos_counter,3) == 0
        fprintf(iges_file,'  %7dP%7d\n',P_num,counter);
        counter=counter+1;
        pos_counter=0;
    end
end

for v_idx=1:length(srf.v_knotvctr)
    pos_counter=pos_counter+1;
    fprintf(iges_file,'%20.12g,',real(srf.v_knotvctr(v_idx)));
    if mod(pos_counter,3) == 0
        fprintf(iges_file,'  %7dP%7d\n',P_num,counter);
        counter=counter+1;
        pos_counter=0;
    end
end

for v_idx=1:srf.v_coef_num
    for u_idx=1:srf.u_coef_num
        pos_counter=pos_counter+1;
        fprintf(iges_file,'%20.12g,',real(weights(u_idx,v_idx,end)));
        if mod(pos_counter,3) == 0
            fprintf(iges_file,'  %7dP%7d\n',P_num,counter);
            counter=counter+1;
            pos_counter=0;
        end
    end
end

for v_idx=1:srf.v_coef_num
    for u_idx=1:srf.u_coef_num
        for dim_idx=1:3
            pos_counter=pos_counter+1;
            fprintf(iges_file,'%20.12g,',real(poles(u_idx,v_idx,dim_idx)));
            if mod(pos_counter,3) == 0
                fprintf(iges_file,'  %7dP%7d\n',P_num,counter);
                counter=counter+1;
                pos_counter=0;
            end
        end
    end
end

% output the ranges
for i=1:4
    pos_counter=pos_counter+1;
    if i == 1
        fprintf(iges_file,'%20.12g,',real(min(srf.u_knotvctr)));
    elseif i == 2
        fprintf(iges_file,'%20.12g,',real(max(srf.u_knotvctr)));
    elseif i == 3
        fprintf(iges_file,'%20.12g,',real(min(srf.v_knotvctr)));
    elseif i == 4
        fprintf(iges_file,'%20.12g;',real(max(srf.v_knotvctr)));
    end
    if mod(pos_counter,3) == 0
        fprintf(iges_file,'  %7dP%7d\n',P_num,counter);
        counter=counter+1;
        pos_counter=0;
    else % we have to close it up anyway
        if i == 4
            for j=1:3-pos_counter
                fprintf(iges_file,'%21s',' ');
            end
            pos_counter=0;
            fprintf(iges_file,'  %7dP%7d\n',P_num,counter);
            counter=counter+1;
        end
    end
end

P_num=P_num+2;
end
