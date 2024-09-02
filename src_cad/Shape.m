classdef Shape < handle
    properties % geometry list
        geom_topo;
        crv_list;
        srf_list;

        crv_name_list;
        srf_name_list;
    end

    methods % define function
        function self=Shape(init_param)
            % generate Shape by IGES file or surface define
            %
            if nargin < 1,init_param=[];end

            if ischar(init_param)
                % generate by IGES file
                iges_filestr=init_param;
                geom_list=readGeomIGES(iges_filestr);

                % load all surface in geom_list
                srf_list=[];srf_name_list={};
                for geom_idx=1:length(geom_list)
                    geom=geom_list{geom_idx};
                    if isa(geom,'Surface')
                        srf_list=[srf_list,geom];
                        srf_name_list=[srf_name_list;{num2str(length(srf_list))}];
                    end
                end

                self.srf_list=srf_list;
                self.srf_name_list=srf_name_list;

                % generate shape topology properties
                %                 self.propagateKnotVectors();
                %                 self.doConnectivity();
            else
                % generate by surface define

            end
        end

        function [crv,crv_idx]=getCurve(self,crv_name)
            % load face from crv_list base on input crv_name
            %
            crv=[];
            crv_idx=[];
            for crv_i=1:length(self.crv_list)
                if strcmpi(self.crv_name_list{crv_i},crv_name)
                    crv=self.srf_list(crv_i);
                    crv_idx=crv_i;
                    return;
                end
            end
        end

        function [srf,srf_idx]=getSurface(self,srf_name)
            % load face from srf_list base on input srf_name
            %
            srf=[];
            srf_idx=[];
            for srf_i=1:length(self.srf_list)
                if strcmpi(self.srf_name_list{srf_i},srf_name)
                    srf=self.srf_list(srf_i);
                    srf_idx=srf_i;
                    return;
                end
            end
        end

        function [grd_list,crd_list]=calGeom(self,param,sym_dir)
            % calculate gird of surface
            %
            if nargin < 3
                sym_dir=[];
                if nargin < 2
                    param=[];
                end
            end
            if isempty(sym_dir),sym_dir=0;end

            % calculate each surface point
            srf_num=length(self.srf_list);
            if sym_dir
                crd_list=repmat(struct('U',[],'V',[]),srf_num*2,1);
                grd_list=repmat(struct('T',[],'X',[],'Y',[],'Z',[],'name',[],'type',[],'orientation',[]),srf_num*2,1);
            else
                crd_list=repmat(struct('U',[],'V',[]),srf_num,1);
                grd_list=repmat(struct('T',[],'X',[],'Y',[],'Z',[],'name',[],'type',[],'orientation',[]),srf_num,1);
            end

            for srf_idx=1:srf_num
                srf=self.srf_list(srf_idx);

                % calculate surface point
                [pnts,us,vs]=srf.calGeom(param);

                grd.T=[];
                grd.X=pnts(:,:,1)';
                grd.Y=pnts(:,:,2)';
                grd.Z=pnts(:,:,3)';
                grd.name=self.srf_name_list{srf_idx};
                grd.type='wgs';
                grd.orientation=0;

                crd.U=us;
                crd.V=vs;

                grd_list(srf_idx)=grd;
                crd_list(srf_idx)=crd;

                % calcualate full scale
                if sym_dir
                    pnts(:,:,sym_dir)=-pnts(:,:,sym_dir);

                    grd.T=[];
                    grd.X=pnts(:,:,1)';
                    grd.Y=pnts(:,:,2)';
                    grd.Z=pnts(:,:,3)';
                    grd.name=self.srf_name_list{srf_idx};
                    grd.type='wgs';
                    grd.orientation=1;

                    crd.U=us;
                    crd.V=vs;

                    grd_list(srf_num+srf_idx)=grd;
                    crd_list(srf_num+srf_idx)=crd;
                end
            end
        end

        function mesh_data=calMesh(self,param,sym_dir)
            % base on calGeom to get LaWGS format mesh data
            %
            if nargin < 3
                sym_dir=[];
                if nargin < 2
                    param=[];
                end
            end

            % generate mesh data
            grd_list=self.calGeom(param,sym_dir);
            srf_num=length(grd_list);
            mesh_data=struct();
            for srf_idx=1:srf_num
                grd=grd_list(srf_idx);
                marker_name=grd.name;

                X=grd.X;
                Y=grd.Y;
                Z=grd.Z;

                % fix normal vector
                if grd.orientation
                    X=flipud(X);
                    Y=flipud(Y);
                    Z=flipud(Z);
                end

                mesh_data.(marker_name).X=X;
                mesh_data.(marker_name).Y=Y;
                mesh_data.(marker_name).Z=Z;
                mesh_data.(marker_name).type='wgs';
            end
        end

        function writeIGES(self,iges_filestr)
            % write shape into IGES file
            %
            if nargin < 2,iges_filestr=[];end
            if isempty(iges_filestr),iges_filestr='geom.iges';end

            writeGeomIGES(num2cell(self.srf_list),iges_filestr);
        end

        function varargout=displayGeom(self,axe_hdl,srf_option,param,sym_dir)
            % draw srf_list on axes handle
            %
            if nargin < 5
                sym_dir=[];
                if nargin < 4
                    param=[];
                    if nargin < 3
                        srf_option=[];
                        if nargin < 2
                            axe_hdl=[];
                        end
                    end
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end
            if isempty(srf_option),srf_option=struct('LineStyle','none','FaceAlpha',0.5);end

            grd_list=self.calGeom(param,sym_dir);

            grd_num=length(grd_list);
            for grd_idx=1:grd_num
                grd=grd_list(grd_idx);

                % draw surface on axe_hdl
                if strcmp(grd.type,'wgs')
                    srf_hdl_list(grd_idx)=surface(axe_hdl,grd.X,grd.Y,grd.Z,srf_option);
                elseif strcmp(grd.type,'stl')
                    srf_hdl_list(grd_idx)=patch(axe_hdl,'faces',trids,'vertices',[grd.X(:),grd.Y(:),grd.Z(:)],'facecolor',grd.Z(:),srf_option);
                else
                    error('unknown grd type');
                end

                % draw name on axe_hdl
                text(axe_hdl,mean(grd.X,"all"),mean(grd.Y,"all"),mean(grd.Z,"all"),grd.name);
            end

            xlabel('x');
            ylabel('y');
            zlabel('z');

            varargout={};
            if nargout > 0
                varargout={srf_hdl_list};
            end
        end

        function varargout=displayDirect(self,axe_hdl)
            % draw direction of surface on axes handle
            %
            if nargin < 2
                axe_hdl=[];
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            srf_num=length(self.srf_list);
            for srf_idx=1:srf_num
                srf=self.srf_list(srf_idx);
                quv_hdl_list{srf_idx}=srf.displayDirect(axe_hdl);
            end
            xlabel('x');
            ylabel('y');
            zlabel('z');

            varargout={};
            if nargout > 0
                varargout={quv_hdl_list};
            end
        end
    end

    methods % topology function
        function self=updateSurfaceCoef(self)
            % Copy the pyGeo list of control points back to the surfaces
            for ii=1:length(self.coef)
                for jj=1:length(self.geom_topo.gIndex{ii})
                    srf_idx=self.geom_topo.gIndex{ii}(jj,1);
                    i=self.geom_topo.gIndex{ii}(jj,2);
                    j=self.geom_topo.gIndex{ii}(jj,3);
                    self.srf_list{srf_idx}.coef(i,j)=double(self.coef(ii));
                end
            end

            for srf_idx=1:length(self.srf_list)
                self.srf_list(srf_idx).setEdgeCurves();
            end
        end

        function self=setSurfaceCoef(self)
            % Set the surface coef list from the pyspline surfaces
            self.coef=zeros(self.geom_topo.nGlobal,3);
            for srf_idx=1:length(self.srf_list)
                srf=self.srf_list(srf_idx);
                for i=1:srf.nCtlu
                    for j=1:srf.nCtlv
                        self.coef(self.topo.lIndex{srf_idx}(i,j),:)=srf.coef(i,j);
                    end
                end
            end
        end

        function doConnectivity(self,fileName,nodeTol,edgeTol)
            % This is the only public edge connectivity function.
            % If fileName exists it loads the file OR it calculates the connectivity
            % and saves to that file.
            %
            % Parameters:
            %   fileName (string): Filename for con file
            %   nodeTol (double): The tolerance for identical nodes
            %   edgeTol (double): The tolerance for midpoint of edges being identical

            srf_num=length(self.srf_list);
            if nargin > 1 && exist(fileName,'file')
                fprintf('Reading Connectivity File: %s\n',fileName);
                self.geom_topo=SurfaceTopology(fileName);
                if ~strcmp(self.initType,'iges')
                    self.propagateKnotVectors();
                end

                sizes=zeros(srf_num,2);
                for isurf=1:srf_num
                    sizes(isurf,:)=[self.srf_list(isurf).nCtlu,self.srf_list(isurf).nCtlv];
                    self.srf_list(isurf).recompute();
                end
                self.geom_topo.calcGlobalNumbering(sizes);
            else
                self.calcConnectivity(nodeTol,edgeTol);
                sizes=zeros(srf_num,2);
                for isurf=1:srf_num
                    sizes(isurf,:)=[self.srf_list(isurf).nCtlu,self.srf_list(isurf).nCtlv];
                end
                self.geom_topo.calcGlobalNumbering(sizes);
                if ~strcmp(self.initType,'iges')
                    self.propagateKnotVectors();
                end
                if nargin > 1
                    fprintf('Writing Connectivity File: %s\n',fileName);
                    self.geom_topo.writeConnectivity(fileName);
                end
            end

            if strcmp(self.initType,'iges')
                self.setSurfaceCoef();
            end
        end

        function calcConnectivity(self,nodeTol,edgeTol)
            % This function attempts to automatically determine the connectivity
            % between the patches

            % Calculate the 4 corners and 4 midpoints for each surface
            srf_num=length(self.srf_list);
            coords=zeros(srf_num,8,3);

            for isurf=1:srf_num
                beg,mid,e=self.srf_list(isurf).getOrigValuesEdge(1);
                coords(isurf,1,:)=beg;
                coords(isurf,2,:)=e;
                coords(isurf,5,:)=mid;
                beg,mid,e=self.srf_list(isurf).getOrigValuesEdge(2);
                coords(isurf,3,:)=beg;
                coords(isurf,4,:)=e;
                coords(isurf,6,:)=mid;
                beg,mid,e=self.srf_list(isurf).getOrigValuesEdge(3);
                coords(isurf,7,:)=mid;
                beg,mid,e=self.srf_list(isurf).getOrigValuesEdge(4);
                coords(isurf,8,:)=mid;
            end

            self.geom_topo=SurfaceTopology(coords,nodeTol,edgeTol);
        end

        function printConnectivity(self)
            % Print the Edge connectivity to the screen
            self.geom_topo.printConnectivity();
        end

        function propagateKnotVectors(self)
            % Propagate the knot vectors to make consistent

            % First get the number of design groups
            nDG=-1;
            ncoef=[];
            for i=1:self.geom_topo.nEdge
                if self.geom_topo.edges(i).dg > nDG
                    nDG=self.geom_topo.edges(i).dg;
                    ncoef(end+1)=self.geom_topo.edges(i).N;
                end
            end
            nDG=nDG+1;

            srf_num=length(self.srf_list);
            for isurf=1:srf_num
                dgU=self.geom_topo.edges(self.geom_topo.edgeLink(isurf,1)).dg;
                dgV=self.geom_topo.edges(self.geom_topo.edgeLink(isurf,3)).dg;
                self.srf_list(isurf).nCtlu=ncoef(dgU);
                self.srf_list(isurf).nCtlv=ncoef(dgV);
                if self.srf_list(isurf).ku < self.srf_list(isurf).nCtlu
                    if self.srf_list(isurf).nCtlu > 4
                        self.srf_list(isurf).ku=4;
                    else
                        self.srf_list(isurf).ku=self.srf_list(isurf).nCtlu;
                    end
                end
                if self.srf_list(isurf).kv < self.srf_list(isurf).nCtlv
                    if self.srf_list(isurf).nCtlv > 4
                        self.srf_list(isurf).kv=4;
                    else
                        self.srf_list(isurf).kv=self.srf_list(isurf).nCtlv;
                    end
                end
                self.srf_list(isurf).calcKnots();
            end

            % Now loop over the number of design groups,accumulate all
            % the knot vectors that corresponds to this dg,then merge them all
            for idg=1:nDG
                knotVectors=[];
                flip=[];
                for isurf=1:srf_num
                    for iedge=1:4
                        if self.geom_topo.edges(self.geom_topo.edgeLink(isurf,iedge)).dg == idg
                            if self.geom_topo.edgeDir(isurf,iedge) == -1
                                flip(end+1)=true;
                            else
                                flip(end+1)=false;
                            end
                            if iedge <= 2
                                knotVec=self.srf_list(isurf).tu;
                            else
                                knotVec=self.srf_list(isurf).tv;
                            end
                            if flip(end)
                                knotVectors(end+1)=fliplr(1 - knotVec);
                            else
                                knotVectors(end+1)=knotVec;
                            end
                        end
                    end
                end

                % Now blend all the knot vectors
                newKnotVec=geo_utils.blendKnotVectors(knotVectors,false);
                newKnotVecFlip=fliplr(1 - newKnotVec);

                counter=1;
                for isurf=1:srf_num
                    for iedge=1:4
                        if self.geom_topo.edges(self.geom_topo.edgeLink(isurf,iedge)).dg == idg
                            if iedge <= 2
                                if flip(counter)
                                    self.srf_list(isurf).tu=newKnotVecFlip;
                                else
                                    self.srf_list(isurf).tu=newKnotVec;
                                end
                            else
                                if flip(counter)
                                    self.srf_list(isurf).tv=newKnotVecFlip;
                                else
                                    self.srf_list(isurf).tv=newKnotVec;
                                end
                            end
                            counter=counter+1;
                        end
                    end
                end
            end
        end

    end

    methods % control function
        function self=translate(self,tran_vctr)
            % translate shape
            %
            for crv_idx=1:length(self.crv_list)
                self.crv_list(crv_idx)=self.crv_list(crv_idx).translate(tran_vctr);
            end

            for srf_idx=1:length(self.srf_list)
                self.srf_list(srf_idx)=self.srf_list(srf_idx).translate(tran_vctr);
            end
        end

        function rotate(self,rot_mat,rot_cntr)
            % rotate shape
            %
            if nargin < 3 || isempty(rot_cntr), rot_cntr=[];end

            for crv_idx=1:length(self.crv_list)
                self.crv_list(crv_idx)=self.crv_list(crv_idx).rotate(rot_mat,rot_cntr);
            end

            for srf_idx=1:length(self.srf_list)
                self.srf_list(srf_idx)=self.srf_list(srf_idx).rotate(rot_mat,rot_cntr);
            end
        end
    end

    methods % parameterization function
        function crd_list=calMeshCoord(self,pnt_list,mkr_idxs_list,sym_dir)
            % calculate all input point local coordinate
            %
            % input:
            % pnt_list(matrix): point_num x dimension matrix
            % mkr_idxs_list(list): struct list, struct is .(marker_name)=pnt_idx_list
            % sym_dir(int): mesh symmetry dimension, 0 is none, 1, 2, 3 is x, y, z
            %
            if nargin < 4
                sym_dir=[];
            end
            if isempty(sym_dir),sym_dir=0;end
            crd_list=struct();

            mkr_name_list=fieldnames(mkr_idxs_list);
            for mkr_idx=1:length(mkr_name_list)
                % get surf
                mkr_name=mkr_name_list{mkr_idx};
                srf=self.getSurface(mkr_name);
                if isempty(srf),continue;end
                pnts_idx=mkr_idxs_list.(mkr_name);

                % calculate coordinate
                if sym_dir
                    % divide mesh into right and left
                    [pnts_idx,pnts_idx_sym,pnts,pnts_sym]=dividePoint(pnts_idx,pnt_list,pnts_sym,0);
                else
                    pnts=pnt_list(pnts_idx,1:3);
                end

                uvs=srf.calCoordinate(pnts);
                crd_list.(mkr_name).index=pnts_idx;
                crd_list.(mkr_name).U=uvs(:,1);
                crd_list.(mkr_name).V=uvs(:,2);

                % calculate symmetry mesh
                if sym_dir
                    pnts_sym(:,sym_dir)=-pnts_sym(:,sym_dir);
                    uvs=srf.calCoordinate(reshape(pnts_sym,size(pnts,1),[],3));
                    sym.index=pnts_idx_sym;
                    sym.U=uvs(:,1);
                    sym.V=uvs(:,2);
                    crd_list.sym=sym;
                end
            end
        end

        function mesh_data=calMeshPoint(self,mkr_crd_list)
            % calculate all mesh point by input coord_list
            %
            % input:
            % coord_list(list): struct list,struct is .(marker_name)=struct(index,U,V)
            %
            mesh_data=struct();

            mkr_name_list=fieldnames(mkr_crd_list);
            for srf_idx=1:length(mkr_name_list)
                % get surface
                mkr_name=mkr_name_list{srf_idx};
                srf=self.getSurface(mkr_name);
                if isempty(srf),continue;end
                pnt_idx=mkr_crd_list.(mkr_name).index;
                U=mkr_crd_list.(mkr_name).U;
                V=mkr_crd_list.(mkr_name).V;

                % calculate coordinate
                pnts=srf.calPoint(U,V);
                mesh_data.(mkr_name).type='scatter';
                mesh_data.(mkr_name).index=pnt_idx;
                mesh_data.(mkr_name).X=pnts(:,:,1);
                mesh_data.(mkr_name).Y=pnts(:,:,2);
                mesh_data.(mkr_name).Z=pnts(:,:,3);
            end
        end
    end
end
