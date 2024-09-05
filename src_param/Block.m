classdef Block < handle
    properties % FFD data
        geom_dim=3; % geometry dimension
        volume_list; % FFD of volume
        volume_name_list; % name of volume
        displace_list; % displacement of control point of FFD
    end

    properties % mesh data
        marker_list; % marker to parameterization list
        mesh_data; % mesh data
        mesh_coord; % mesh loacl coordination in FFD
    end

    methods % define function
        function self=Block(init_param)
            % generate Block
            %
            if isa(init_param,'Volume') % initialize FFD by input volume list
                vol_list=init_param(:);
                vol_num=length(vol_list);
                vol_name_list=strsplit(num2str(1:vol_num,'B%d '))';
            elseif isa(init_param,'char') % initialize FFD by input FFD define file
                FFD_filestr=init_param;

                % check if input file exist
                if ~exist(FFD_filestr,'file')
                    error('Block: cfg file do not exist');
                end

                vol_list=Block.readFFD(FFD_filestr);
                vol_num=length(vol_list);
                vol_name_list=strsplit(num2str(1:vol_num,'B%d '))';
            end

            % initialize displacement list
            disp_list=cell(vol_num,1);
            for vol_idx=1:vol_num
                vol=vol_list(vol_idx);
                pole_num=numel(vol.coefs);
                disp=zeros(pole_num,1);

                disp_list{vol_idx}=disp;
            end

            self.volume_list=vol_list;
            self.volume_name_list=vol_name_list;
            self.displace_list=disp_list;
        end

        function setPoleDisplace(self,vol_idxs,is,js,ks,ds,vals)
            % modify FFD control point displacement of Block
            %
            vol_list=self.volume_list;
            disp_list=self.displace_list;
            vol_num=length(vol_list);

            for vol_idx=1:vol_num
                vol=vol_list(vol_idx);
                disp=disp_list{vol_idx};
                idx=find(vol_idxs == vol_idx);
                sub=sub2ind(size(vol.coefs),is(idx),js(idx),ks(idx),ds(idx));
                disp(sub)=vals(idx);

                disp_list{vol_idx}=disp;
            end

            self.displace_list=disp_list;
        end

        function writeBlock(self,FFD_filestr)
            % write block define into file
            %
            if nargin < 2
                FFD_filestr=[];
            end

            vol_list=self.volume_list;
            Block.writeFFD(vol_list,FFD_filestr)
        end

        function displayBlock(self,axe_hdl,pole_option)
            % display all FFD on axes handle
            %
            if nargin < 3
                pole_option=[];
                if nargin < 2
                    axe_hdl=[];
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            vol_list=self.volume_list;
            vol_name_list=self.volume_name_list;
            disp_list=self.displace_list;
            vol_num=length(vol_list);

            for vol_idx=1:vol_num
                vol=vol_list(vol_idx);
                disp=disp_list{vol_idx};
                vol_name=vol_name_list{vol_idx};

                vol.coefs=vol.coefs+reshape(disp,size(vol.coefs)); % add displacement of control point

                vol.displayPole(axe_hdl,pole_option);
                vol.displayDirect(axe_hdl);

                poles=vol.getPoles();
                cntr_pnt=[mean(poles(:,:,:,1),'all'),mean(poles(:,:,:,2),'all'),mean(poles(:,:,:,3),'all')];
                text(cntr_pnt(1),cntr_pnt(2),cntr_pnt(3),vol_name);
            end

            axis equal;view(3);

            x_range=xlim();
            y_range=ylim();
            z_range=zlim();
            center=[mean(x_range),mean(y_range),mean(z_range)];
            range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
            xlim([center(1)-range,center(1)+range]);
            ylim([center(2)-range,center(2)+range]);
            zlim([center(3)-range,center(3)+range]);
        end
    end

    methods % parameterization function
        function mesh_coord=attachMesh(self,mesh_data,marker_list)
            % calculate local parameter of point on FFD
            %
            % input:
            % mesh_data(struct list): data about mesh
            % marker_list(cell): mesh marker to parameterization
            %
            if nargin < 3
                marker_list=[];
            end
            if isempty(marker_list),marker_list=fieldnames(mesh_data);end
            for mkr_idx=1:length(marker_list),if strcmp(marker_list{mkr_idx},'geometry'),marker_list(mkr_idx)=[];break;end;end % remove geometry

            self.mesh_data=mesh_data;

            mesh_coord=struct();
            for mkr_idx=1:length(marker_list)
                mkr_name=marker_list{mkr_idx};
                mkr=mesh_data.(mkr_name);

                if strcmp(mkr.type,'wgs')

                elseif strcmp(mkr.type,'stl')
                    pnt_list=mkr.element_list;

                    [vol_idx_list,uvw_list,~,in_vol_bools,~]=self.attachPoint(pnt_list);
                    index=find(in_vol_bools);

                    mesh_coord.(mkr_name).FFD_ID=vol_idx_list(index);
                    mesh_coord.(mkr_name).index=index;
                    mesh_coord.(mkr_name).U=uvw_list(index,1);
                    mesh_coord.(mkr_name).V=uvw_list(index,2);
                    mesh_coord.(mkr_name).W=uvw_list(index,3);
                elseif strcmp(mkr.type,'mesh')

                end
            end

            self.marker_list=marker_list;
            self.mesh_coord=mesh_coord;
        end

        function [mesh_data,marker_list]=updataMesh(self)
            % calculate new point coordination by FFD
            %
            marker_list=self.marker_list;
            mesh_data=self.mesh_data;

            for mkr_idx=1:length(marker_list)
                mkr_name=marker_list{mkr_idx};
                mkr=mesh_data.(mkr_name);
                crd=self.mesh_coord.(mkr_name);

                if strcmp(mkr.type,'wgs')

                elseif strcmp(mkr.type,'stl')
                    pnt_list=mkr.element_list;

                    vol_idx_list=crd.FFD_ID;
                    uvw_list=[crd.U,crd.V,crd.W];
                    index=crd.index;

                    pnts=self.updataPoint(vol_idx_list,uvw_list);
                    pnt_list(index,:)=pnts;
                    mkr.element_list=pnt_list;
                elseif strcmp(mkr.type,'mesh')

                end

                mesh_data.(mkr_name)=mkr;
            end
        end

        function [vol_idx_list,uvw_list,dis_list,in_vol_bools,in_cvh_bools]=attachPoint(self,pnts)
            % attach point to volume list
            %
            vol_list=self.volume_list;
            vol_num=length(vol_list);
            atch_tol=sqrt(eps);

            pnt_num=size(pnts,1);
            vol_idx_list=zeros(pnt_num,1,'int32');
            uvw_list=nan(pnt_num,3);
            dis_list=1e10*ones(pnt_num,3);
            in_vol_bools=false(pnt_num,1);
            in_cvh_bools=false(pnt_num,1);

            for vol_idx=1:vol_num
                vol=vol_list(vol_idx);
                poles=vol.getPoles();
                poles=reshape(poles,[],3);

                % first step, check if inside convex hull of all volume
                hull=convhull(poles,'Simplify',true);

                % calculate hyperplane of convex hull
                vctr_i=poles(hull(:,2),:)-poles(hull(:,1),:);
                vctr_j=poles(hull(:,3),:)-poles(hull(:,2),:);
                nmvctr=cross(vctr_i,vctr_j,2);
                offset=-dot(nmvctr,poles(hull(:,1),:),2);
                dis_plane=pnts*nmvctr'+offset';

                % point in convex hull of this volume
                in_vol=all(dis_plane <= eps,2);
                in_conv_idx=find(in_vol);

                in_cvh_bools=in_cvh_bools | in_vol;

                % second step, project to volume and check if inside volume
                uvw_vol=vol.projectPoint(pnts(in_conv_idx,:),1e-10);
                pnts_proj=vol.calPoint(uvw_vol);
                dis_vol=sum(abs(pnts(in_conv_idx,:)-pnts_proj),2);
                dis_list(in_conv_idx)=dis_vol;

                % point in volume
                in_dis=dis_vol < atch_tol;
                in_vol(in_conv_idx)=in_vol(in_conv_idx) & in_dis;

                vol_idx_list(in_vol,:)=vol_idx;
                uvw_list(in_vol,:)=uvw_vol(in_dis,:);
                in_vol_bools=in_vol_bools | in_vol;
            end
        end

        function pnts=updataPoint(self,vol_idx_list,uvw_list)
            % calculate new point coordination by local parameter
            %
            vol_list=self.volume_list;
            disp_list=self.displace_list;
            vol_num=length(vol_list);

            pnt_num=size(uvw_list,1);
            pnts=zeros(pnt_num,self.geom_dim);

            for vol_idx=1:vol_num
                vol=vol_list(vol_idx);
                disp=disp_list{vol_idx};
                vol.coefs=vol.coefs+reshape(disp,size(vol.coefs)); % add displacement of control point

                idx=find(vol_idx_list == vol_idx);
                pnts(idx,:)=vol.calPoint(uvw_list(idx,:));
            end
        end
    end

    methods(Static) % IO function
        function vol_list=readFFD(FFD_filestr)
            % read FFD define from file
            %
            FFD_file=fopen(FFD_filestr,'r');

            vol_num=textscan(FFD_file,'%d\n',1);vol_num=vol_num{1};
            vol_shape=textscan(FFD_file,'%d %d %d',vol_num);vol_shape=[vol_shape{:}];
            vol_list=Volume.empty(vol_num,0);
            for vol_idx=1:vol_num
                pole_num=vol_shape(vol_idx,1)*vol_shape(vol_idx,2)*vol_shape(vol_idx,3);

                X=textscan(FFD_file,'%f',pole_num);X=X{1};
                Y=textscan(FFD_file,'%f',pole_num);Y=Y{1};
                Z=textscan(FFD_file,'%f',pole_num);Z=Z{1};

                X=reshape(X,vol_shape(vol_idx,:));
                Y=reshape(Y,vol_shape(vol_idx,:));
                Z=reshape(Z,vol_shape(vol_idx,:));

                poles=cat(4,X,Y,Z);

                u_deg=min(3,size(poles,1)-1);
                v_deg=min(3,size(poles,2)-1);
                w_deg=min(3,size(poles,3)-1);

                vol_list(vol_idx)=Volume(poles,u_deg,v_deg,w_deg); % create volume list by poles
            end
            fclose(FFD_file);
            clear('FFD_file');
        end

        function writeFFD(vol_list,FFD_filestr)
            % write FFD poles define into file
            %
            if nargin < 2
                FFD_filestr=[];
            end
            if isempty(FFD_filestr),FFD_filestr='FFDxyz.dat';end

            FFD_file=fopen(FFD_filestr,'w');
            vol_num=length(vol_list);
            fprintf(FFD_file,'%d\n',vol_num);

            poles_list=cell(vol_num,1);
            for vol_idx=1:vol_num
                vol=vol_list(vol_idx);
                poles=vol.getPoles();
                FFD_shape=size(poles);
                fprintf(FFD_file,'%d %d %d\n',FFD_shape(1:end-1));
                poles_list{vol_idx}=poles;
            end

            for vol_idx=1:vol_num
                poles=poles_list{vol_idx};
                out_format=[repmat('%.15e ',1,size(poles,1)-1),'%.15e\n'];
                fprintf(FFD_file,out_format,poles);
            end

            fclose(FFD_file);
            clear('FFD_file');
        end
    end
end
