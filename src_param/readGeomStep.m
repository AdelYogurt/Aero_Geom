function geom_list=readGeomStep(step_filestr)
% read topology entity from step file
%

step_file=fopen(step_filestr,'r');
str_all=fread(step_file)';
str_all=char(str_all);
fclose(step_file);
clear('step_file');

% extract data
idx=strfind(str_all,'DATA;');
str_all=str_all(idx+5:end);
idx=strfind(str_all,'ENDSEC;');
str_all=str_all(1:idx-1);
str_all(str_all == ' ')='';
str_all(str_all == newline)='';
str_all(str_all == 13)='';
Data=textscan(str_all,'#%d=%s','Delimiter',';');

% sort data for index
data_num=max(Data{1});
Didx=(1:data_num)';
Dstr=repmat({''},data_num,1);
Didx(Data{1})=Data{1};
Dstr(Data{1})=Data{2};
clear('Data');clear('str_all');

data_num=length(Didx);
% find out MANIFOLD_SOLID_BREP or SHELL_BASED_SURFACE_MODEL to start
Shl_didx=[];
for str_idx=1:data_num
    str=Dstr{str_idx};

    % check and solid shell
    if contains(str,'SHELL_BASED_SURFACE_MODEL')
        idx=find(str == ',',1,"first");
        str=str(idx+2:end-2);
        shl_didx=textscan(str,'#%d','Delimiter',',');
        shl_didx=cell2mat(shl_didx);
        Shl_didx=[Shl_didx;shl_didx];
        continue;
    end

    if contains(str,'MANIFOLD_SOLID_BREP')
        idx=find(str == ',',1,"first");
        str=str(idx+1:end-1);
        shl_didx=textscan(str,'#%d','Delimiter',',');
        shl_didx=cell2mat(shl_didx);
        Shl_didx=[Shl_didx;shl_didx];
        continue;
    end
end

% read CLOSED/OPEN SHELL
shl_num=length(Shl_didx);
geom_list=repmat(Shell('',[]),shl_num,1);
for shl_idx=1:shl_num
    shl_didx=Shl_didx(shl_idx);
    shl=readShell(Didx,Dstr,shl_didx);
    geom_list(shl_idx)=shl;
end

end

%% topology step

function shl=readShell(Didx,Dstr,shl_didx)
% generate shell entity
%
str=Dstr{shl_didx};
% process shell string, CLOSED/OPEN SHELL

% read shell name face
idx=find(str == '(',1,"first");
str=str(idx+1:end-1);
idx=find(str == ',',1,"first");
shl_name=str(1:idx-1);

% read face index
str=str(idx+1:end);
Fce_didx=textscan(str(2:end-1),'#%d','Delimiter',',');
Fce_didx=cell2mat(Fce_didx);

% read ADVANCED_FACE
fce_num=length(Fce_didx);
fce_list=repmat(Face('',[]),fce_num,1);
for fce_idx=1:fce_num
    fce_didx=Fce_didx(fce_idx);
    fce=readFace(Didx,Dstr,fce_didx);
    fce_list(fce_idx)=fce;
end

shl=Shell(shl_name,fce_list);
end

function fce=readFace(Didx,Dstr,fce_didx)
% generate face entity
%
str=Dstr{fce_didx};
% process face string, ADVANCED_FACE

% read face name
idx=find(str == '(',1,"first");
str=str(idx+1:end-1);
idx=find(str == ',',1,"first");
fce_name=str(1:idx-1);
str=str(idx+1:end);

% read FACE OUTER/ BOUND index
idx=find(str == ')',1,"first");
Wir_didx=textscan(str(2:idx-1),'#%d','Delimiter',',');
Wir_didx=cell2mat(Wir_didx);
str=str(idx+2:end);

% read SURFACE index
idx=find(str == ',',1,"first");
srf_didx=textscan(str(1:idx-1),'#%d');
srf_didx=cell2mat(srf_didx);

% read FACE OUTER/ BOUND
wir_num=length(Wir_didx);
wir_list=repmat(Wire('',[]),wir_num,1);
for wir_idx=1:wir_num
    wir_didx=Wir_didx(wir_idx);
    wir=readLoop(Didx,Dstr,wir_didx);
    wir_list(wir_idx)=wir;
end

% read SURFACE
srf=readSurface(Didx,Dstr,srf_didx,wir_list);

fce=Face(fce_name,srf,wir_list);
end

function [wir,otr]=readLoop(Didx,Dstr,wir_didx)
% generate wire entity
%
str=Dstr{wir_didx};

% FACE_OUTER_BOUND/FACE_BOUND
if contains(str,'FACE_OUTER_BOUND'), otr=true;
elseif contains(str,'FACE_BOUND'), otr=false;
else, error('readLoop: unknown BOUND type');
end

idx=find(str == ',',1,"first");
str=str(idx+1:end-1);

% read EDGE_LOOP index
idx=find(str == ',',1,"first");
wir_didx=textscan(str(1:idx-1),'#%d','Delimiter',',');
wir_didx=cell2mat(wir_didx);
str=str(idx+1:end);

% read reverse boolean
if contains(str,'F'), wir_reverse=true;
else, wir_reverse=false; end
    
% EDGE_LOOP
str=Dstr{wir_didx};

idx=find(str == '(',1,"first");
str=str(idx+1:end-1);

% read wire name
idx=find(str == ',',1,"first");
wir_name=str(1:idx-1);
str=str(idx+1:end);

% read ORIENTED_EDGE index
idx=find(str == ')',1,"first");
Edg_didx=textscan(str(2:idx-1),'#%d','Delimiter',',');
Edg_didx=cell2mat(Edg_didx);

edg_num=length(Edg_didx);
edg_list=repmat(Edge('',[]),edg_num,1);
% ORIENTED_EDGE
for edg_idx=1:edg_num
    idx_edge=Edg_didx(edg_idx);
    edge=readEdge(Didx,Dstr,idx_edge);
    edg_list(edg_idx)=edge;
end

wir=Wire(wir_name,edg_list);
if wir_reverse, wir.reverse();end
end

function edg=readEdge(Didx,Dstr,edg_didx)
% generate edge entity
%
str=Dstr{edg_didx};

idx=find(str == ',');
str=str(idx(3)+1:end-1);

% read EDGE_CURVE index
idx=find(str == ',',1,"first");
edg_didx=textscan(str(1:idx-1),'#%d','Delimiter',',');
edg_didx=cell2mat(edg_didx);
str=str(idx+1:end);

% read reverse boolean
if contains(str,'F'), edg_reverse=true;
else, edg_reverse=false; end

% EDGE_LOOP
str=Dstr{edg_didx};

% read edge name
idx=find(str == '(',1,"first");
str=str(idx+1:end-1);
idx=find(str == ',',1,"first");
edg_name=str(1:idx-1);
str=str(idx+1:end);

% read VERTEX index
idx=find(str == ',');
Vtx_didx=textscan(str(1:idx(2)-1),'#%d','Delimiter',',');
Vtx_didx=cell2mat(Vtx_didx);
str=str(idx(2)+1:end);

% read CURVE index
idx=find(str == ',',1,"first");
crv_didx=textscan(str(1:idx-1),'#%d','Delimiter',',');
crv_didx=cell2mat(crv_didx);

% read VERTEX
vtx_list=readVertex(Didx,Dstr,Vtx_didx);

crv=readCurve(Didx,Dstr,crv_didx,vtx_list);

edg=Edge(edg_name,crv,vtx_list);
if edg_reverse, edg.reverse();end
end

function Vtx=readVertex(Didx,Dstr,Vtx_didx)
% generate vertex entity
%
str=Dstr(Vtx_didx);
str=strjoin(str,'\n');
Pnt_didx=textscan(str,'VERTEX_POINT(%s#%d)','Delimiter',{',','\n'});
Pnt_didx=Pnt_didx{2};
Vtx=readPoint(Didx,Dstr,Pnt_didx);
end

%% geometry step

function srf=readSurface(Didx,Dstr,srf_didx,wir_list)
% read surface
%
str=Dstr{srf_didx};

if str(1) == '('
    % if is list
    str=str(2:end-1);
    if contains(str,'RATIONAL_B_SPLINE_SURFACE')
        % scan str to find str data
        while ~isempty(str)
            sidx=find(str == '(',1,'first');
            eidx=matchParenthesis(str,sidx);
            if strcmp(str(1:sidx-1),'B_SPLINE_SURFACE')
                str_pole=str(sidx+1:eidx-1);
            elseif strcmp(str(1:sidx-1),'B_SPLINE_SURFACE_WITH_KNOTS')
                str_knot=str(sidx+1:eidx-1);
            elseif strcmp(str(1:sidx-1),'RATIONAL_B_SPLINE_SURFACE')
                str_weight=str(sidx+1:eidx-1);
            end
            str=str(eidx+1:end);
        end

        % read UDegree, VDegree
        idx=find(str_pole == ',',2,"first");
        Degree=textscan(str_pole,'%f%f','Delimiter',',');
        UDegree=Degree{1};VDegree=Degree{2};
        str_pole=str_pole(idx(2)+1:end);

        % read Poles index
        idx=matchParenthesis(str_pole,1);
        data_pole=str_pole(2:(idx-2));
        data_pole=strsplit(data_pole,'),');
        Pnt_didx=[];
        for u_idx=1:length(data_pole)
            pnt_didx=textscan(data_pole{u_idx}(2:end),'#%d','Delimiter',',');
            Pnt_didx=[Pnt_didx,cell2mat(pnt_didx)];
        end
        str_pole=str_pole(idx+1:end);

        % read curve_form, u_closed, v_closed, self_intersect
        srf_data=textscan(str_pole,'%s%s%s%s','Delimiter',',');
        surface_form=srf_data{1};
        u_closed=srf_data{2};
        if contains(u_closed,'T'), u_closed=true;
        elseif contains(u_closed,'F'), u_closed=false;end
        v_closed=srf_data{3};
        if contains(v_closed,'T'), v_closed=true;
        elseif contains(v_closed,'F'), v_closed=false;end
        self_intersect=srf_data{4};

        % read Poles
        [v_num,u_num]=size(Pnt_didx);
        Poles=readPoint(Didx,Dstr,Pnt_didx(:));
        Poles=reshape(Poles,v_num,u_num,[]);

        % read UMults
        idx=find(str_knot == ')',1,"first")+1;
        UMults=textscan(str_knot(2:(idx-2)),'%f','Delimiter',',');
        UMults=cell2mat(UMults);
        str_knot=str_knot(idx+1:end);

        % read VMults
        idx=find(str_knot == ')',1,"first")+1;
        VMults=textscan(str_knot(2:(idx-2)),'%f','Delimiter',',');
        VMults=cell2mat(VMults);
        str_knot=str_knot(idx+1:end);

        % read UKnots
        idx=find(str_knot == ')',1,"first")+1;
        UKnots=textscan(str_knot(2:(idx-2)),'%f','Delimiter',',');
        UKnots=cell2mat(UKnots);UKnots=(UKnots-min(UKnots))/(max(UKnots)-min(UKnots));
        str_knot=str_knot(idx+1:end);

        % read VKnots
        idx=find(str_knot == ')',1,"first")+1;
        VKnots=textscan(str_knot(2:(idx-2)),'%f','Delimiter',',');
        VKnots=cell2mat(VKnots);VKnots=(VKnots-min(VKnots))/(max(VKnots)-min(VKnots));

        % read Weights
        idx=matchParenthesis(str_weight,1);
        data_weight=str_weight(2:(idx-2));
        data_weight=strsplit(data_weight,'),');
        Weights=[];
        for u_idx=1:length(data_weight)
            weight=textscan(data_weight{u_idx}(2:end),'%f','Delimiter',',');
            Weights=[Weights,cell2mat(weight)];
        end

        knot_spec=[];
    else
        error('readSurface: unsupport surface type');
    end
else
    if contains(str,'B_SPLINE_SURFACE_WITH_KNOTS')
        % read UDegree, VDegree
        idx=find(str == ',',3,"first");
        Degree=textscan(str((idx(1)+1):(idx(3)-1)),'%f%f','Delimiter',',');
        UDegree=Degree{1};VDegree=Degree{2};
        str=str(idx(3)+1:end);

        % read Poles index
        idx=matchParenthesis(str,1);
        data_pole=str(2:(idx-2));
        data_pole=strsplit(data_pole,'),');
        Pnt_didx=[];
        for u_idx=1:length(data_pole)
            pnt_didx=textscan(data_pole{u_idx}(2:end),'#%d','Delimiter',',');
            Pnt_didx=[Pnt_didx,cell2mat(pnt_didx)];
        end
        str=str(idx+2:end);

        % read curve_form, u_closed, v_closed, self_intersect
        idx=find(str == ',',4,"first");
        srf_data=textscan(str(1:(idx(4)-1)),'%s%s%s%s','Delimiter',',');
        surface_form=srf_data{1};
        u_closed=srf_data{2};
        if contains(u_closed,'T'), u_closed=true;
        elseif contains(u_closed,'F'), u_closed=false;end
        v_closed=srf_data{3};
        if contains(v_closed,'T'), v_closed=true;
        elseif contains(v_closed,'F'), v_closed=false;end
        self_intersect=srf_data{4};
        str=str(idx(4)+1:end);

        % read UMults
        idx=find(str == ')',1,"first")+1;
        UMults=textscan(str(2:(idx-2)),'%f','Delimiter',',');
        UMults=cell2mat(UMults);
        str=str(idx+1:end);

        % read VMults
        idx=find(str == ')',1,"first")+1;
        VMults=textscan(str(2:(idx-2)),'%f','Delimiter',',');
        VMults=cell2mat(VMults);
        str=str(idx+1:end);

        % read UKnots
        idx=find(str == ')',1,"first")+1;
        UKnots=textscan(str(2:(idx-2)),'%f','Delimiter',',');
        UKnots=cell2mat(UKnots);UKnots=(UKnots-min(UKnots))/(max(UKnots)-min(UKnots));
        str=str(idx+1:end);

        % read VKnots
        idx=find(str == ')',1,"first")+1;
        VKnots=textscan(str(2:(idx-2)),'%f','Delimiter',',');
        VKnots=cell2mat(VKnots);VKnots=(VKnots-min(VKnots))/(max(VKnots)-min(VKnots));
        str=str(idx+1:end);

        % read knot_spec
        knot_spec=str(1:end-1);

        % read Poles
        [v_num,u_num]=size(Pnt_didx);
        Poles=readPoint(Didx,Dstr,Pnt_didx(:));
        Poles=reshape(Poles,v_num,u_num,[]);
        
        Weights=[];
    elseif contains(str,'PLANE')
        edg_list=wir_list(1).edge_list;

        switch length(edg_list)
            case 2
                crv_u0=edg_list(1).curve;
                crv_1v=[];
                crv_u1=copy(edg_list(2).curve).reverse();
                crv_0v=[];
            case 3
                crv_u0=edg_list(1).curve;
                crv_1v=edg_list(2).curve;
                crv_u1=copy(edg_list(3).curve).reverse();
                crv_0v=[];
            case 4
                crv_u0=edg_list(1).curve;
                crv_1v=edg_list(2).curve;
                crv_u1=copy(edg_list(3).curve).reverse();
                crv_0v=copy(edg_list(4).curve).reverse();
            otherwise
                idx=find(str == ',',1,"first");
                str=str(idx+1:end);

                % read axis index and radius
                ci_data=textscan(str(1:end-1),'#%d%f','Delimiter',',');
                ax_didx=ci_data{1};radius=ci_data{2};
                ax=readAxis(Didx,Dstr,ax_didx);
                pnt_cen=ax(1,:);
                zdir=ax(2,:)/norm(ax(2,:));
                xdir=ax(3,:)/norm(ax(3,:));
                ydir=cross(zdir,xdir);

                % findout bound of wir, and create plane directly
                bou_box=wir_list(1).bound_box;
                pln_bou=(bou_box-pnt_cen)*[xdir;ydir;zdir]';
                pln_bou=[min(pln_bou,[],1);max(pln_bou,[],1)];
                Poles=[pln_bou(1,1),pln_bou(1,2),0;
                    pln_bou(1,1),pln_bou(2,2),0;
                    pln_bou(2,1),pln_bou(1,2),0;
                    pln_bou(2,1),pln_bou(2,2),0];

                Poles=Poles*[xdir;ydir;zdir]+pnt_cen;
                Poles=reshape(Poles,2,2,[]);
                srf=Surface(Poles);
                return;
        end

        srf=boundCurveToSurface(crv_u0,crv_u1,crv_0v,crv_1v);
        return;
    elseif contains(str,'SPHERICAL_SURFACE') || contains(str,'CYLINDRICAL_SURFACE')
        edg_list=wir_list(1).edge_list;

        % if edge num equal to 4, edge order is 0v, u1, 1v, u0
        switch length(edg_list)
            case 2
                crv_u0=edg_list(1).curve;
                crv_1v=[];
                crv_u1=copy(edg_list(2).curve).reverse();
                crv_0v=[];
            case 3
                crv_u0=edg_list(1).curve;
                crv_1v=edg_list(2).curve;
                crv_u1=copy(edg_list(3).curve).reverse();
                crv_0v=[];
            case 4
                crv_u0=edg_list(1).curve;
                crv_1v=edg_list(2).curve;
                crv_u1=copy(edg_list(3).curve).reverse();
                crv_0v=copy(edg_list(4).curve).reverse();
        end

        srf=boundCurveToSurface(crv_u0,crv_u1,crv_0v,crv_1v);
        return;
    else
        error('readSurface: unsupport surface type');
    end
end

srf=Surface(Poles,UDegree,VDegree,UMults,VMults,UKnots,VKnots,Weights);

srf.surface_form=surface_form;
srf.UPeriodic=u_closed;
srf.VPeriodic=v_closed;
srf.self_intersect=self_intersect;
srf.knot_spec=knot_spec;
end

function crv=readCurve(Didx,Dstr,crv_didx,vtx_list)
% read curve
%
str=Dstr{crv_didx};

if str(1) == '('
    % if is list
    if contains(str,'RATIONAL_B_SPLINE_CURVE')
        % read B_SPLINE_SURFACE
        sidx=strfind(str,'B_SPLINE_CURVE');
        sidx=sidx(str(sidx+length('B_SPLINE_CURVE')) == '(' &...
            str(sidx-1) ~= '_')+length('B_SPLINE_CURVE');
        eidx=matchParenthesis(str,sidx);
        str_pole=str(sidx+1:eidx-1);

        % read B_SPLINE_SURFACE_WITH_KNOTS
        sidx=strfind(str,'B_SPLINE_CURVE_WITH_KNOTS')+length('B_SPLINE_CURVE_WITH_KNOTS');
        eidx=matchParenthesis(str,sidx);
        str_knot=str(sidx+1:eidx-1);

        % read RATIONAL_B_SPLINE_CURVE
        sidx=strfind(str,'RATIONAL_B_SPLINE_CURVE')+length('RATIONAL_B_SPLINE_CURVE');
        eidx=matchParenthesis(str,sidx);
        str_weight=str(sidx+1:eidx-1);

        % read Degree
        idx=find(str_pole == ',',1,"first");
        Degree=textscan(str_pole(1:(idx(1)-1)),'%f');
        Degree=cell2mat(Degree);
        str_pole=str_pole(idx(1)+1:end);

        % read Poles index
        idx=find(str_pole == ')',1,"first")+1;
        Pnt_didx=textscan(str_pole(2:(idx-2)),'#%d','Delimiter',',');
        Pnt_didx=cell2mat(Pnt_didx);
        str_pole=str_pole(idx+1:end);

        % read Poles
        Poles=readPoint(Didx,Dstr,Pnt_didx);

        % read curve_form, closed_curve, self_intersect
        crv_data=textscan(str_pole,'%s%s%s','Delimiter',',');
        curve_form=crv_data{1};
        closed_curve=crv_data{2};
        if contains(closed_curve,'T'), closed_curve=true;
        elseif contains(closed_curve,'F'), closed_curve=false;end
        self_intersect=crv_data{3};

        % read Mults
        idx=find(str_knot == ')',1,"first")+1;
        Mults=textscan(str_knot(2:(idx-2)),'%f','Delimiter',',');
        Mults=cell2mat(Mults);
        str_knot=str_knot(idx+1:end);

        % read Knots
        idx=find(str_knot == ')',1,"first")+1;
        Knots=textscan(str_knot(2:(idx-2)),'%f','Delimiter',',');
        Knots=cell2mat(Knots);Knots=(Knots-min(Knots))/(max(Knots)-min(Knots));

        % read Weights
        idx=find(str_weight == ')',1,"first")+1;
        Weights=textscan(str_weight(2:(idx-2)),'%f','Delimiter',',');
        Weights=cell2mat(Weights);

        knot_spec=[];
    else
        error('readCurve: unsupport curve type');
    end
else
    if contains(str,'B_SPLINE_CURVE')
        % read Degree
        idx=find(str == ',',2,"first");
        Degree=textscan(str((idx(1)+1):(idx(2)-1)),'%f');
        Degree=cell2mat(Degree);
        str=str(idx(2)+1:end);

        % read Poles index
        idx=find(str == ')',1,"first")+1;
        Pnt_didx=textscan(str(2:(idx-2)),'#%d','Delimiter',',');
        Pnt_didx=cell2mat(Pnt_didx);
        str=str(idx+1:end);

        % read curve_form, closed_curve, self_intersect
        idx=find(str == ',',3,"first");
        crv_data=textscan(str(1:(idx(3)-1)),'%s%s%s','Delimiter',',');
        curve_form=crv_data{1};
        closed_curve=crv_data{2};
        if contains(closed_curve,'T'), closed_curve=true;
        elseif contains(closed_curve,'F'), closed_curve=false;end
        self_intersect=crv_data{3};
        str=str(idx(3)+1:end);

        % read Mults
        idx=find(str == ')',1,"first")+1;
        Mults=textscan(str(2:(idx-2)),'%f','Delimiter',',');
        Mults=cell2mat(Mults);
        str=str(idx+1:end);

        % read Knots
        idx=find(str == ')',1,"first")+1;
        Knots=textscan(str(2:(idx-2)),'%f','Delimiter',',');
        Knots=cell2mat(Knots);Knots=(Knots-min(Knots))/(max(Knots)-min(Knots));
        str=str(idx+1:end);

        % read knot_spec
        knot_spec=str(1:end-1);

        % read Poles
        Poles=readPoint(Didx,Dstr,Pnt_didx);

        Weights=[];
    elseif contains(str,'LINE')
        crv=Curve(vtx_list);
        return;
    elseif contains(str,'CIRCLE')
        idx=find(str == ',',1,"first");
        str=str(idx+1:end);

        % read axis index and radius
        ci_data=textscan(str(1:end-1),'#%d%f','Delimiter',',');
        ax_didx=ci_data{1};radius=ci_data{2};

        % read axis
        ax=readAxis(Didx,Dstr,ax_didx);
        pnt_cen=ax(1,:);
        zdir=ax(2,:)/norm(ax(2,:));
        xdir=ax(3,:)/norm(ax(3,:));
        ydir=cross(zdir,xdir);

        svtr=vtx_list(1,:)-pnt_cen;svtr=svtr/norm(svtr);
        evtr=vtx_list(2,:)-pnt_cen;evtr=evtr/norm(evtr);
        svtr=svtr*[xdir;ydir;zdir]';
        evtr=evtr*[xdir;ydir;zdir]';
        sang=cart2pol(svtr(1),svtr(2));eang=cart2pol(evtr(1),evtr(2));

        crv=circle(radius,sang,eang,[pnt_cen;xdir;ydir]);
        return;
    else
        error('readCurve: unsupport curve type');
    end
end

crv=Curve(Poles,Degree,Mults,Knots,Weights);
crv.curve_form=curve_form;
crv.Periodic=closed_curve;
crv.self_intersect=self_intersect;
crv.knot_spec=knot_spec;
end

function ax=readAxis(Didx,Dstr,ax_didx)
% read axis
%
str=Dstr{ax_didx};

% read Poles index
idx=find(str == ',',1,"first");
Pnt_didx=textscan(str(idx+1:end-1),'#%d','Delimiter',',');
Pnt_didx=cell2mat(Pnt_didx);

ax=readPoint(Didx,Dstr,Pnt_didx);
end

function Pnt=readPoint(Didx,Dstr,Pnt_didx)
% read point
%
str=Dstr(Pnt_didx);
str=strjoin(str,'\n');
str(str == '(' | str == ')')='';
Pnt=textscan(str,'%s%f%f%f','Delimiter',{',','\n'});
Pnt=cell2mat(Pnt(:,2:4));
end

%% unit

function eidx=matchParenthesis(str,eidx)
% str(eidx) must be (
% 
stack=1;
while eidx < length(str) && stack
    eidx=eidx+1;
    if str(eidx) == '('
        stack=stack+1;
    end
    if str(eidx) == ')'
        stack=stack-1;
    end
end
end