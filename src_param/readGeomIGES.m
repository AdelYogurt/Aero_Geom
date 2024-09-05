function geom_list=readGeomIGES(iges_filestr)
% read geometry entity from IGES file
%
iges_file=fopen(iges_filestr,'r');
str_all=fread(iges_file)';
str_all=char(str_all);
fclose(iges_file);
clear('iges_file');

str_all=strrep(str_all,';',',');str_all=strrep(str_all,char(13),'');
datstr=strsplit(str_all,'\n');datstr=datstr(1:end-1);

srt_num=str2double(datstr{end}(2:8));
gnrl_num=str2double(datstr{end}(10:16));
dir_num=str2double(datstr{end}(18:24));
param_num=str2double(datstr{end}(26:32));

% now we know how many lines we have to deal with
dir_offset=srt_num+gnrl_num;
para_offset=dir_offset+dir_num;

% read geometry head
datidx=zeros(dir_num/2,3);
enty_num=0;
idxlin=[];idxrevl=[];idxcrv=[];idxsrf=[];
% directory lines is a multiple of 2
for i=1:2:dir_num
    type_idx=str2double(datstr{i+dir_offset}(2:8));
    start=str2double(datstr{i+dir_offset}(10:16));
    num_lines=str2double(datstr{i+1+dir_offset}(26:32));
    datidx(enty_num+1,:)=[type_idx,start,num_lines];

    switch type_idx
%         case 110 % 110 is line
%             idxlin=[idxlin,geom_num+1];
        case 120
            idxrevl=[idxrevl,enty_num+1];
%         case 126 % 126 is spline curve
%             idxcrv=[idxcrv,geom_num+1];
        case 128 % 128 is spline surface
            idxsrf=[idxsrf,enty_num+1];
    end
    enty_num=enty_num+1;
end
fprintf("Found %d geometry entity in IGES File.\n",enty_num);

% read geometry parameter
geom_num=length(idxlin)+length(idxrevl)+length(idxcrv)+length(idxsrf);
geom_list=cell(1,geom_num);
param_offset=dir_offset+dir_num;

geom_read=0;
for idx=1:length(idxlin)
    lin_idx=idxlin(idx);
    lin=readLineParam(datidx,datstr,lin_idx,param_offset);
    geom_list{geom_read+1}=lin;
    geom_read=geom_read+1;
end

for idx=1:length(idxrevl)
    revl_idx=idxrevl(idx);
    srf=readRevolveParam(datidx,datstr,revl_idx,param_offset);
    geom_list{geom_read+1}=srf;
    geom_read=geom_read+1;
end

for idx=1:length(idxcrv)
    crv_idx=idxcrv(idx);
    crv=readCurveParam(datidx,datstr,crv_idx,param_offset);
    geom_list{geom_read+1}=crv;
    geom_read=geom_read+1;
end

for idx=1:length(idxsrf)
    srf_idx=idxsrf(idx);
    srf=readSurfaceParam(datidx,datstr,srf_idx,param_offset);
    geom_list{geom_read+1}=srf;
    geom_read=geom_read+1;
end

fprintf("Successful load %d geometry entity in IGES File.\n",geom_num);
end

%% geometry IGES

function srf=readSurfaceParam(datidx,datstr,srf_idx,param_offset)
% read surface parameter information from IGES file
%
offset=datidx(srf_idx,2)+param_offset-1;

% read all string data of surface
str_num=datidx(srf_idx,3);
str_srf=[datstr{(1:datidx(srf_idx,3))+offset}];
idx=(1:65)'+(0:(str_num-1))*80;
str_srf=str_srf(idx(:));

% convert string into numeric
data=textscan(str_srf,' %f,');
data=data{1};

% extract surface data
n_number=round(data(2))+1;
v_coef_num=round(data(3))+1;
u_order=round(data(4));
v_order=round(data(5));

u_closed=data(6);
v_closed=data(7);

counter=11;
u_knotvctr=data(counter:counter+n_number+u_order)';
counter=counter+n_number+u_order+1;

v_knotvctr=data(counter:counter+v_coef_num+v_order)';
counter=counter+v_coef_num+v_order+1;

weights=data(counter:counter+n_number*v_coef_num-1);
weights=reshape(weights,n_number,v_coef_num);
counter=counter+n_number*v_coef_num;

poles=data(counter:counter+n_number*v_coef_num*3-1);
poles=reshape(poles,3,n_number,v_coef_num);poles=permute(poles,[2,3,1]);
counter=counter+n_number*v_coef_num*3;

coefs=cat(3,poles,weights);

srf=Surface(coefs,{u_knotvctr,v_knotvctr});
end

function crv=readCurveParam(datidx,datstr,crv_idx,param_offset)
% read curve parameter information from IGES file
%
offset=datidx(crv_idx,2)+param_offset-1;

% read all string data of surface
str_num=datidx(crv_idx,3);
str_crv=[datstr{(1:datidx(crv_idx,3))+offset}];
idx=(1:65)'+(0:(str_num-1))*80;
str_crv=str_crv(idx(:));

% convert string into numeric
data=textscan(str_crv,' %f,');
data=data{1};

% extract curve data
u_coef_num=round(data(2))+1;
u_order=round(data(3));

planed=data(4);
u_closed=data(5);

counter=8;
u_knotvctr=data(counter:counter+u_coef_num+u_order)';
counter=counter+u_coef_num+u_order+1;

weights=data(counter:counter+u_coef_num-1);
weights=reshape(weights,u_coef_num,1);
counter=counter+u_coef_num;

poles=data(counter:counter+u_coef_num*3-1);
poles=reshape(poles,3,u_coef_num);poles=permute(poles,[2,1]);
counter=counter+u_coef_num*3;

coefs=cat(2,poles,weights);

crv=Curve(coefs,u_knotvctr);
end

function srf=readRevolveParam(datidx,datstr,revl_idx,param_offset)
% read curve parameter information from IGES file
%
offset=datidx(revl_idx,2)+param_offset-1;

% read all string data of surface
str_num=datidx(revl_idx,3);
str_revl=[datstr{(1:datidx(revl_idx,3))+offset}];
idx=(1:65)'+(0:(str_num-1))*80;
str_revl=str_revl(idx(:));

% convert string into numeric
data=textscan(str_revl,' %f,');
data=data{1};

axs_idx=(data(2)-1)/2+1;
geom_idx=(data(3)-1)/2+1;
srt_ang=data(4);
eng_ang=data(5);

% rotation coordination
axs=readLineParam(datidx,datstr,axs_idx,param_offset);
axs=axs.getPoles();
pnt_rot=axs(1,:);
vctr_rot=axs(2,:)-axs(1,:);vctr_rot=vctr_rot/norm(vctr_rot);

% rotation geometry
type_idx=datidx(geom_idx,1);
switch type_idx
    case 110 % 110 is line
        crv=readLineParam(datidx,datstr,geom_idx,param_offset);
    case 126 % 126 is spline curve
        crv=readCurveParam(datidx,datstr,geom_idx,param_offset);
end

crv=crv.rotate(axang2rotm([vctr_rot,srt_ang]),pnt_rot);
srf=constrSurfaceRevolve(crv,pnt_rot,vctr_rot,srt_ang,eng_ang);
end

function lin=readLineParam(datidx,datstr,lin_idx,param_offset)
% read curve parameter information from IGES file
%
offset=datidx(lin_idx,2)+param_offset-1;

% read all string data of surface
str_num=datidx(lin_idx,3);
str_lin=[datstr{(1:datidx(lin_idx,3))+offset}];
idx=(1:65)'+(0:(str_num-1))*80;
str_lin=str_lin(idx(:));

% convert string into numeric
data=textscan(str_lin,' %f,');
data=data{1};

% extract line data
poles=data(2:7);
poles=reshape(poles,3,2);
poles=poles';

lin=Curve(poles);
end
