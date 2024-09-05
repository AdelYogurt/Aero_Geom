classdef Volume
    % Non-Uniform Rational B-Splines Volume
    % define reference to step standard
    %
    properties % Explicit Attributes
        form='.B-NURBS.'; % surface_form
        coef_dim=[]; % dimension of coefs
        u_coef_num; % u number of coefs
        v_coef_num; % v number of coefs
        w_coef_num; % w number of coefs
        coefs; % coefs matrix, equal to permute(cat(4,poles.*weights,weights),[4,3,2,1]);
        u_knotvctr; % u_knot_vector
        v_knotvctr; % v_knot_vector
        w_knotvctr; % v_knot_vector
        u_order=[]; % u_degree
        v_order=[]; % v_degree
        w_order=[]; % w_degree
    end

    properties % Derivate Attributes
        % step standard properties
        knot_spec=''; % knot_spec
        u_closed=false; % UPeriodic (boolean)
        v_closed=false; % VPeriodic (boolean)
        w_closed=false; % WPeriodic (boolean)
        intersected=false; % self_intersect (boolean)

        u_deriv_vol=Volume.empty(0);
        v_deriv_vol=Volume.empty(0);
        w_deriv_vol=Volume.empty(0);
    end

    methods % define volume
        function self=Volume(poles,u_degree,v_degree,w_degree,u_mults,v_mults,w_mults,u_knots,v_knots,w_knots,weights)
            % generate Non-Uniform Rational B-Splines Volume
            %
            % calling:
            % vol=Volume(poles,u_degree,v_degree,w_degree,u_mults,v_mults,w_mults,u_knots,v_knots,w_knots,weights)
            % vol=Volume(coefs,{u_knotvctr,v_knotvctr,w_knotvctr})
            %
            % input:
            % poles (matrix): control point, u_pole_num x v_pole_num x w_pole_num x dimension matrix
            % u_degree (matrix): optional input
            % v_degree (matrix): optional input
            % w_degree (matrix): optional input
            % u_mults (matrix): optional input
            % v_mults (matrix): optional input
            % w_mults (matrix): optional input
            % u_knots (matrix): optional input
            % v_knots (matrix): optional input
            % w_knots (matrix): optional input
            % weights (matrix): optional input
            %
            % output:
            % Volume
            %
            % notice:
            % degree default is pole_num-1, which will be Bezier volume
            %
            if nargin < 11
                weights=[];
                if nargin < 10
                    w_knots=[];
                    if nargin < 9
                        v_knots=[];
                        if nargin < 8
                            u_knots=[];
                            if nargin < 7
                                w_mults=[];
                                if nargin < 6
                                    v_mults=[];
                                    if nargin < 5
                                        u_mults=[];
                                        if nargin < 4
                                            w_degree=[];
                                            if nargin < 3
                                                v_degree=[];
                                                if nargin < 2
                                                    u_degree=[];
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            [u_pole_num,v_pole_num,w_pole_num,coef_dim]=size(poles);

            % default value
            if isempty(u_degree),u_degree=u_pole_num-1;end
            if isempty(v_degree),v_degree=v_pole_num-1;end
            if isempty(w_degree),w_degree=w_pole_num-1;end

            % different input mode
            if iscell(u_degree) && length(u_degree) > 1
                coefs=poles;
                [u_knotvrtc,v_knotvrtc,w_knotvrtc]=u_degree{:};
                u_knotvrtc=sort(u_knotvrtc);
                v_knotvrtc=sort(v_knotvrtc);
                w_knotvrtc=sort(w_knotvrtc);
                u_order=length(u_knotvrtc)-u_pole_num-1;
                v_order=length(v_knotvrtc)-v_pole_num-1;
                w_order=length(w_knotvrtc)-w_pole_num-1;

                if u_pole_num < (u_order+1) || v_pole_num < (v_order+1) || w_pole_num < (w_order+1)
                    error('Volume: pole_num less than order+1');
                end
            elseif length(u_degree) == 1
                % u default value
                if isempty(u_mults) && isempty(u_knots)
                    u_mults=[u_degree+1,ones(1,u_pole_num-u_degree-1),u_degree+1];
                    u_knots=linspace(0,1,u_pole_num-u_degree+1);
                elseif ~isempty(u_mults) && isempty(u_knots)
                    u_knots=linspace(0,1,length(u_mults));
                elseif isempty(u_mults) && ~isempty(u_knots)
                    error('Volume: need u_mults input');
                end
                u_knotvrtc=baseKnotVctr(u_mults,u_knots);
                u_knotvrtc=sort(u_knotvrtc);

                % v default value
                if isempty(v_mults) && isempty(v_knots)
                    v_mults=[v_degree+1,ones(1,v_pole_num-v_degree-1),v_degree+1];
                    v_knots=linspace(0,1,v_pole_num-v_degree+1);
                elseif ~isempty(mults) && isempty(v_knots)
                    v_knots=linspace(0,1,length(v_mults));
                elseif isempty(v_mults) && ~isempty(v_knots)
                    error('Volume: need v_mults input');
                end
                v_knotvrtc=baseKnotVctr(v_mults,v_knots);
                v_knotvrtc=sort(v_knotvrtc);

                % w default value
                if isempty(w_mults) && isempty(w_knots)
                    w_mults=[w_degree+1,ones(1,w_pole_num-w_degree-1),w_degree+1];
                    w_knots=linspace(0,1,w_pole_num-w_degree+1);
                elseif ~isempty(mults) && isempty(w_knots)
                    w_knots=linspace(0,1,length(w_mults));
                elseif isempty(w_mults) && ~isempty(w_knots)
                    error('Volume: need w_mults input');
                end
                w_knotvrtc=baseKnotVctr(w_mults,w_knots);
                w_knotvrtc=sort(w_knotvrtc);

                if isempty(weights), weights=ones(u_pole_num,v_pole_num,w_pole_num,1);end
                if ~isempty(weights), weights=reshape(weights,[u_pole_num,v_pole_num,w_pole_num,1]);end

                if u_pole_num < (u_degree+1) || v_pole_num < (v_degree+1) || w_pole_num < (w_degree+1)
                    error('Volume: pole_num less than degree+1');
                end

                if length(u_knotvrtc) ~= u_pole_num+u_degree+1 || length(v_knotvrtc) ~= v_pole_num+v_degree+1 || length(w_knotvrtc) ~= w_pole_num+w_degree+1
                    error('Volume: knot_num is not equal to pole_num+degree+1');
                end

                u_order=u_degree;
                v_order=v_degree;
                w_order=w_degree;
                coefs=cat(4,poles.*weights,weights);
                coef_dim=coef_dim+1;
            else
                error('Volume: error input format');
            end

            self.coef_dim=coef_dim;
            self.u_coef_num=u_pole_num;
            self.v_coef_num=v_pole_num;
            self.w_coef_num=w_pole_num;
            self.u_order=u_order;
            self.v_order=v_order;
            self.w_order=w_order;
            self.coefs=coefs;
            self.u_knotvctr=u_knotvrtc;
            self.v_knotvctr=v_knotvrtc;
            self.w_knotvctr=w_knotvrtc;
        end

        function self=deriv(self,deriv_time)
            % generate derivate BSpline for calculate gradient
            %

            % recursion to derivate data
            if deriv_time > 0
                if isempty(self.u_deriv_vol)
                    % generate derivate volume along the u direction
                    coefs_mat=permute(self.coefs,[1,2,3,4]);
                    coefs_mat=reshape(coefs_mat,[self.u_coef_num,self.v_coef_num*self.w_coef_num*self.coef_dim]);
                    [coefs_mat,~,u_knotvrtc]=BSpline.deriv(coefs_mat,self.u_order,self.u_knotvctr);
                    coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.v_coef_num,self.w_coef_num,self.coef_dim]);
                    coefs_mat=permute(coefs_mat,[1,2,3,4]);

                    self.u_deriv_vol=Volume(coefs_mat,{u_knotvrtc,self.v_knotvctr,self.w_knotvctr});
                end

                if isempty(self.v_deriv_vol)
                    % generate derivate volume along the v direction
                    coefs_mat=permute(self.coefs,[2,1,3,4]);
                    coefs_mat=reshape(coefs_mat,[self.v_coef_num,self.u_coef_num*self.w_coef_num*self.coef_dim]);
                    [coefs_mat,~,v_knotvrtc]=BSpline.deriv(coefs_mat,self.v_order,self.v_knotvctr);
                    coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.u_coef_num,self.w_coef_num,self.coef_dim]);
                    coefs_mat=permute(coefs_mat,[2,1,3,4]);

                    self.v_deriv_vol=Volume(coefs_mat,{self.u_knotvctr,v_knotvrtc,self.w_knotvctr});
                end

                if isempty(self.w_deriv_vol)
                    % generate derivate volume along the w direction
                    coefs_mat=permute(self.coefs,[3,2,1,4]);
                    coefs_mat=reshape(coefs_mat,[self.w_coef_num,self.v_coef_num*self.u_coef_num*self.coef_dim]);
                    [coefs_mat,~,w_knotvrtc]=BSpline().deriv(coefs_mat,self.w_order,self.w_knotvctr);
                    coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.v_coef_num,self.u_coef_num,self.coef_dim]);
                    coefs_mat=permute(coefs_mat,[3,2,1,4]);

                    self.w_deriv_vol=Volume(coefs_mat,{self.u_knotvctr,self.v_knotvctr,w_knotvrtc});
                end

                % donot output data to avoid matlab repeat deriv_crv
                self.v_deriv_vol=self.v_deriv_vol.deriv(deriv_time-1);
                self.u_deriv_vol=self.u_deriv_vol.deriv(deriv_time-1);
                self.w_deriv_vol=self.w_deriv_vol.deriv(deriv_time-1);
            end
        end

        function [poles,weights]=getPoles(self)
            % convert coefs to poles and weight
            %
            poles=self.coefs(:,:,:,1:end-1);
            weights=self.coefs(:,:,:,end);
            weights(weights == 0)=1;
            poles=poles./weights;
        end

        function [pnts,wts]=calPoint(self,us,vs,ws)
            % calculate point on volume
            %
            % calling:
            % [pnts,wts]=vol.calPoint({u_list,v_list,w_list})
            % [pnts,wts]=vol.calPoint(uvw_list)
            % [pnts,wts]=srf.calPoint(us,vs,ws)
            %
            % notice:
            % pnts and wts is ndgrid generate format, [u,v,w]
            % last dimension of uvws should equal to 3
            %
            if nargin < 3,vs=[];ws=[];end
            if isempty(vs),uvws=us;end

            if iscell(us)
                us=uvws{1}(:);vs=uvws{2}(:);ws=uvws{3}(:);
            else
                uvws=cat(4,us,vs,ws);
                siz=size(uvws);uvws=reshape(uvws,[],3);
                us=uvws(:,1);vs=uvws(:,2);ws=uvws(:,3);
            end

            [u_N_list,u_idx_srt,~]=baseFcnN(us,self.u_order,self.u_knotvctr);
            [v_N_list,v_idx_srt,~]=baseFcnN(vs,self.v_order,self.v_knotvctr);
            [w_N_list,w_idx_srt,~]=baseFcnN(ws,self.w_order,self.w_knotvctr);

            if iscell(uvws)
                % input is {u_list,v_list,w_list}, mesh point
                u_num=length(us);v_num=length(vs);w_num=length(ws);

                % evaluate along the u direction
                coefs_u=reshape(self.coefs,[self.u_coef_num,self.v_coef_num*self.w_coef_num*self.coef_dim]);
                pnts=zeros(length(us),size(coefs_u,2));
                for deg_idx=1:self.u_order+1
                    pnts=pnts+u_N_list(:,deg_idx).*coefs_u(u_idx_srt+(deg_idx-1),:);
                end
                pnts=reshape(pnts,[u_num,self.v_coef_num,self.w_coef_num,self.coef_dim]);

                % evaluate along the v direction
                coefs_v=permute(pnts,[2,1,3,4]);
                coefs_v=reshape(coefs_v,[self.v_coef_num,self.coef_dim*u_num*self.w_coef_num]);
                pnts=zeros(length(vs),size(coefs_v,2));
                for deg_idx=1:self.v_order+1
                    pnts=pnts+v_N_list(:,deg_idx).*coefs_v(v_idx_srt+(deg_idx-1),:);
                end
                pnts=reshape(pnts,[v_num,u_num,self.w_coef_num,self.coef_dim]);
                pnts=permute(pnts,[2,1,3,4]);

                % evaluate along the w direction
                coefs_w=permute(pnts,[3,2,1,4]);
                coefs_w=reshape(coefs_w,[self.w_coef_num,v_num*u_num*self.coef_dim]);
                pnts=zeros(length(ws),size(coefs_w,2));
                for deg_idx=1:self.w_order+1
                    pnts=pnts+w_N_list(:,deg_idx).*coefs_w(w_idx_srt+(deg_idx-1),:);
                end
                pnts=reshape(pnts,[w_num,v_num,u_num,self.coef_dim]);
                pnts=permute(pnts,[3,2,1,4]);

                wts=pnts(:,:,:,end);
                pnts=pnts(:,:,:,1:end-1);
                if nargout < 2
                    pnts=pnts./wts;
                end
            else
                % input is uvw_list, scatter point
                pnt_num=size(uvws,1);

                % evaluate along the u direction
                coefs_u=reshape(self.coefs,[self.u_coef_num,self.v_coef_num*self.w_coef_num*self.coef_dim]);
                pnts_u=zeros(length(us),size(coefs_u,2));
                for deg_idx=1:self.u_order+1
                    pnts_u=pnts_u+u_N_list(:,deg_idx).*coefs_u(u_idx_srt+(deg_idx-1),:);
                end

                % evaluate along the v direction
                coefs_v=reshape(pnts_u,[],self.w_coef_num*self.coef_dim);
                pnts_v=zeros(pnt_num,self.w_coef_num*self.coef_dim);
                for deg_idx=1:self.v_order+1
                    idx=sub2ind([pnt_num,self.v_coef_num],(1:pnt_num)',v_idx_srt+(deg_idx-1));
                    pnts_v=pnts_v+v_N_list(:,deg_idx).*coefs_v(idx,:);
                end

                % evaluate along the w direction
                coefs_w=reshape(pnts_v,[],self.coef_dim);
                pnts=zeros(pnt_num,self.coef_dim);
                for deg_idx=1:self.w_order+1
                    idx=sub2ind([pnt_num,self.w_coef_num],(1:pnt_num)',w_idx_srt+(deg_idx-1));
                    pnts=pnts+w_N_list(:,deg_idx).*coefs_w(idx,:);
                end

                wts=pnts(:,end);
                pnts=pnts(:,1:end-1);
                if nargout < 2
                    pnts=pnts./wts;
                end

                pnts=reshape(pnts,[siz(1:end-1),self.coef_dim-1]);
                wts=reshape(wts,[siz(1:end-1),1]);
            end
        end

        function varargout=calGradient(self,us,vs,ws)
            % calculate gradient of volume
            %
            % calling:
            % vol.calGradient({u_list,v_list,w_list})
            % vol.calGradient(uvw_list)
            %
            if nargin < 4,ws=[];if nargin < 3,vs=[];end;end
            deriv_time=nargout-1;
            self=self.deriv(deriv_time);

            % recursion to calculate gradient
            if deriv_time > 0
                % calculate initial curve point and weight
                [c0,w0]=self.calPoint(us,vs,ws);
                p0=c0./w0;

                % calculate first derivative
                [c1u,w1u]=self.u_deriv_vol.calPoint(us,vs,ws);
                [c1v,w1v]=self.v_deriv_vol.calPoint(us,vs,ws);
                [c1w,w1w]=self.w_deriv_vol.calPoint(us,vs,ws);
                p1u=(c1u-p0.*w1u)./w0;
                p1v=(c1v-p0.*w1v)./w0;
                p1w=(c1w-p0.*w1w)./w0;
                p1={p1u,p1v,p1w};

                varargout={p0,p1};
                if nargout > 2
                    % calculate second derivative
                    [c2uu,w2uu]=self.u_deriv_vol.u_deriv_vol.calPoint(us,vs,ws);
                    [c2uv,w2uv]=self.u_deriv_vol.v_deriv_vol.calPoint(us,vs,ws);
                    [c2uw,w2uw]=self.u_deriv_vol.w_deriv_vol.calPoint(us,vs,ws);
                    [c2vv,w2vv]=self.v_deriv_vol.v_deriv_vol.calPoint(us,vs,ws);
                    [c2vw,w2vw]=self.v_deriv_vol.w_deriv_vol.calPoint(us,vs,ws);
                    [c2ww,w2ww]=self.w_deriv_vol.w_deriv_vol.calPoint(us,vs,ws);

                    p2{1,1}=(c2uu-2*p1{1}.*w1u)./w0-p0.*w2uu;
                    p2{1,2}=(c2uv-p1u.*w1v-p1v.*w1u-p0.*w2uv)./w0;
                    p2{1,3}=(c2uw-p1u.*w1w-p1w.*w1u-p0.*w2uw)./w0;
                    p2{2,1}=p2{1,2};
                    p2{2,2}=(c2vv-2*p1{2}.*w1v)./w0-p0.*w2vv;
                    p2{2,3}=(c2vw-p1v.*w1w-p1w.*w1v-p0.*w2vw)./w0;
                    p2{3,1}=p2{1,3};
                    p2{3,2}=p2{2,3};
                    p2{3,3}=(c2ww-2*p1{3}.*w1w)./w0-p0.*w2ww;

                    varargout=[varargout,{p2}];
                end
            end
        end

        function varargout=displayGeom(self,axe_hdl,vol_option,param)
            % draw volume on axes handle
            %
            if nargin < 4
                param=[];
                if nargin < 3
                    vol_option=[];
                    if nargin < 2
                        axe_hdl=[];
                    end
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % default draw option
            if isempty(vol_option), vol_option=struct('LineStyle','none','FaceAlpha',0.5);end

            % draw Points on axe_hdl
            srf_list=self.extract();
            for srf_idx=1:length(srf_list)
                vol_hdl(srf_idx)=srf_list(srf_idx).displayGeom(axe_hdl,vol_option,param);
            end
            xlabel('\itX');
            ylabel('\itY');
            zlabel('\itZ');
            
            if nargout > 0,varargout={vol_hdl};end
        end

        function varargout=displayPole(self,axe_hdl,option)
            % draw curve on figure handle
            %
            if nargin < 3
                option=[];
                if nargin < 2
                    axe_hdl=[];
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % default draw option
            surface_option=struct('Marker','s','MarkerEdgeColor','r','EdgeColor','r','LineStyle','--','FaceAlpha',0);
            if ~isempty(option)
                names=fieldnames(option);
                for idx=1:length(names),surface_option.(names{idx})=option.(names{idx});end
            end

            % draw poles on axe_hdl
            [poles,~]=self.getPoles();
            for k_idx=1:self.w_coef_num
                vol_hdl(k_idx)=surface(axe_hdl,squeeze(poles(:,:,k_idx,1))',squeeze(poles(:,:,k_idx,2))',squeeze(poles(:,:,k_idx,3))',surface_option);
            end
            for i_idx=1:self.u_coef_num
                vol_hdl(i_idx+self.w_coef_num)=surface(axe_hdl,squeeze(poles(i_idx,:,:,1))',squeeze(poles(i_idx,:,:,2))',squeeze(poles(i_idx,:,:,3))',surface_option);
            end
            xlabel('\itX');
            ylabel('\itY');
            zlabel('\itZ');
            
            if nargout > 0,varargout={vol_hdl};end
        end

        function varargout=displayDirect(self,axe_hdl)
            % draw curve on figure handle
            %
            if nargin < 2
                axe_hdl=[];
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % calculate coord on curve
            pnt_mesh=self.coefs(1:2,1:2,1:2,1:end-1);
            wt_mesh=self.coefs(1:2,1:2,1:2,end);
            wt_mesh(wt_mesh == 0)=1;
            pnt_mesh=pnt_mesh./wt_mesh;
            origin=pnt_mesh(1,1,1,:);
            coord=squeeze([pnt_mesh(2,1,1,:)-origin,pnt_mesh(1,2,1,:)-origin,pnt_mesh(1,1,2,:)-origin]);
            coord=(coord./vecnorm(coord,2,2))';

            % scale of quiver
            scale=norm(diff([axe_hdl.XLim;axe_hdl.YLim;axe_hdl.ZLim],1,2));
            coord=coord*scale*0.05;

            % draw coord on axe_hdl
            if self.coef_dim-1 == 3
                hold on;
                quv_hdl(1)=quiver3(axe_hdl,origin(1),origin(2),origin(3),coord(1,1),coord(2,1),coord(3,1),0,'Color','r','LineWidth',1);
                quv_hdl(2)=quiver3(axe_hdl,origin(1),origin(2),origin(3),coord(1,2),coord(2,2),coord(3,2),0,'Color','g','LineWidth',1);
                quv_hdl(3)=quiver3(axe_hdl,origin(1),origin(2),origin(3),coord(1,3),coord(2,3),coord(3,3),0,'Color','b','LineWidth',1);
                hold off;
                xlabel('\itX');
                ylabel('\itY');
                zlabel('\itZ');
            end

            if nargout > 0,varargout={quv_hdl};end
        end

        function varargout=extract(self,sides)
            % get boundary surface
            %
            if nargin < 2,sides=[];end
            if isempty(sides), sides=1:6;end

            srf_list=Surface.empty(0,6);
            for idx=1:3
                srf_idx1=(idx-1)*2+1;
                srf_idx2=(idx-1)*2+2;

                switch idx
                    case 1 % U
                        coefs1=squeeze(self.coefs(1,:,:,:));
                        coefs2=squeeze(self.coefs(end,:,:,:));
                        srf_list(srf_idx1)=Surface(coefs1,{self.v_knotvctr,self.w_knotvctr});
                        srf_list(srf_idx2)=Surface(coefs2,{self.v_knotvctr,self.w_knotvctr});
                    case 2 % V
                        coefs1=squeeze(self.coefs(:,1,:,:));
                        coefs2=squeeze(self.coefs(:,end,:,:));
                        srf_list(srf_idx1)=Surface(coefs1,{self.u_knotvctr,self.w_knotvctr});
                        srf_list(srf_idx2)=Surface(coefs2,{self.u_knotvctr,self.w_knotvctr});
                    case 3 % W
                        coefs1=squeeze(self.coefs(:,:,1,:));
                        coefs2=squeeze(self.coefs(:,:,end,:));
                        srf_list(srf_idx1)=Surface(coefs1,{self.u_knotvctr,self.v_knotvctr});
                        srf_list(srf_idx2)=Surface(coefs2,{self.u_knotvctr,self.v_knotvctr});
                end
            end

            srf_list=srf_list(sides);
            if nargout > 1, varargout=num2cell(srf_list);
            else, varargout={srf_list};end
        end
    end

    methods % control volume
        function self=permute(self,ord)
            % change volume parameter order
            %
            knots_list={self.u_knotvctr,self.v_knotvctr,self.w_knotvctr};
            self=Volume(permute(self.coefs,[ord,4]),knots_list(ord));
        end

        function self=reverseU(self)
            % revese U direction of volume
            %
            self.coefs=self.coefs(end:-1:1,:,:,:);
            min_k=min(self.u_knotvctr);max_k=max(self.u_knotvctr);dk=max_k-min_k;
            self.u_knotvctr=dk-(self.u_knotvctr(end:-1:1)-min_k)+min_k;
        end

        function self=reverseV(self)
            % revese V direction of volume
            %
            self.coefs=self.coefs(:,end:-1:1,:,:);
            min_k=min(self.v_knotvctr);max_k=max(self.v_knotvctr);dk=max_k-min_k;
            self.v_knotvctr=dk-(self.v_knotvctr(end:-1:1)-min_k)+min_k;
        end

        function self=reverseW(self)
            % revese W direction of volume
            %
            self.coefs=self.coefs(:,:,end:-1:1,:);
            min_k=min(self.w_knotvctr);max_k=max(self.w_knotvctr);dk=max_k-min_k;
            self.w_knotvctr=dk-(self.w_knotvctr(end:-1:1)-min_k)+min_k;
        end

        function self=addOrder(self,u_times,v_times,w_times)
            % increase volume degree
            %
            if u_times <= 0 && v_times <= 0 && w_times <= 0, return;end

            % modify order along the u direction
            if u_times <= 0
                coefs_mat=permute(self.coefs,[1,2,3,4]);
                coefs_mat=reshape(coefs_mat,[self.u_coef_num,self.v_coef_num*self.w_coef_num*self.coef_dim]);
                [coefs_mat,self.u_order,self.u_knotvctr]=BSpline.addOrder(coefs_mat,self.u_order,self.u_knotvctr,u_times);
                coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.v_coef_num,self.w_coef_num,self.coef_dim]);
                self.coefs=permute(coefs_mat,[1,2,3,4]);
                self.w_coef_num=size(self.coefs,4);
            end

            % modify order along the v direction
            if v_times <= 0
                coefs_mat=permute(self.coefs,[2,1,3,4]);
                coefs_mat=reshape(coefs_mat,[self.v_coef_num,self.u_coef_num*self.w_coef_num*self.coef_dim]);
                [coefs_mat,self.v_order,self.v_knotvctr]=BSpline.addOrder(coefs_mat,self.v_order,self.v_knotvctr,v_times);
                coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.u_coef_num,self.w_coef_num,self.coef_dim]);
                self.coefs=permute(coefs_mat,[2,1,3,4]);
                self.v_coef_num=size(self.coefs,3);
            end

            % modify order along the w direction
            if w_times <= 0
                coefs_mat=permute(self.coefs,[3,2,1,4]);
                coefs_mat=reshape(coefs_mat,[self.w_coef_num,self.v_coef_num*self.u_coef_num*self.coef_dim]);
                [coefs_mat,self.w_order,self.w_knotvctr]=BSpline.addOrder(coefs_mat,self.w_order,self.w_knotvctr,w_times);
                coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.v_coef_num,self.u_coef_num,self.coef_dim]);
                self.coefs=permute(coefs_mat,[3,2,1,4]);
                self.u_coef_num=size(self.coefs,2);
            end
        end

        function self=insertKnot(self,u_knot_ins,v_knot_ins,w_knot_ins)
            % insert knot to volume
            %
            if isempty(u_knot_ins) && isempty(v_knot_ins) && isempty(w_knot_ins), return;end

            % modify order along the u direction
            if ~isempty(u_knot_ins)
                coefs_mat=permute(self.coefs,[1,2,3,4]);
                coefs_mat=reshape(coefs_mat,[self.u_coef_num,self.v_coef_num*self.w_coef_num*self.coef_dim]);
                [coefs_mat,self.u_order,self.u_knotvctr]=BSpline.insertKnot(coefs_mat,self.u_order,self.u_knotvctr,u_knot_ins);
                coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.v_coef_num,self.w_coef_num,self.coef_dim]);
                self.coefs=permute(coefs_mat,[1,2,3,4]);
                self.w_coef_num=size(self.coefs,4);
            end

            % insert knots along the v direction
            if ~isempty(v_knot_ins)
                coefs_mat=permute(self.coefs,[1,2,4,3]);
                coefs_mat=reshape(coefs_mat,[self.v_coef_num,self.u_coef_num*self.w_coef_num*self.coef_dim]);
                [coefs_mat,self.v_order,self.v_knotvctr]=BSpline.insertKnot(coefs_mat,self.v_order,self.v_knotvctr,v_knot_ins);
                coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.u_coef_num,self.w_coef_num,self.coef_dim]);
                self.coefs=permute(coefs_mat,[1,2,4,3]);
                self.v_coef_num=size(self.coefs,3);
            end

            % insert knots along the w direction
            if ~isempty(w_knot_ins)
                coefs_mat=permute(self.coefs,[1,3,4,2]);
                coefs_mat=reshape(coefs_mat,[self.w_coef_num,self.v_coef_num*self.u_coef_num*self.coef_dim]);
                [coefs_mat,self.w_order,self.w_knotvctr]=BSpline.insertKnot(coefs_mat,self.w_order,self.w_knotvctr,w_knot_ins);
                coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.v_coef_num,self.u_coef_num,self.coef_dim]);
                self.coefs=permute(coefs_mat,[1,4,2,3]);
                self.u_coef_num=size(self.coefs,2);
            end
        end

        function self=translate(self,tran_vctr)
            % translate volume
            %
            tran_vctr=reshape(tran_vctr,1,1,1,[]);
            self.coefs(:,:,:,1:end-1)=self.coefs(:,:,:,1:end-1)+self.coefs(:,:,:,end).*tran_vctr(1:self.coef_dim-1);
        end

        function self=rotate(self,rot_mat,rot_cntr)
            % rotate volume
            %
            if nargin < 3 || isempty(rot_cntr), rot_cntr=zeros(self.coef_dim-1,1);end
            self=self.translate(-rot_cntr);
            self.coefs(:,:,:,1:end-1)=reshape(reshape(self.coefs(:,:,:,1:end-1),[],self.coef_dim-1)*rot_mat',self.u_coef_num,self.v_coef_num,self.w_coef_num,self.coef_dim-1);
            self=self.translate(rot_cntr);
        end
    end

    methods % calculate coord
        function uvw_list=calCoordinate(self,pnts_init,geom_tol)
            % base on X, Y, Z calculate local coordinate in volume
            %
            if nargin < 3, geom_tol=[];end
            if isempty(geom_tol), geom_tol=sqrt(eps);end

            % find point to start
            uvw_list=self.findNearest(pnts_init,20);

            % use project function to adjust parameter
            uvw_list=self.projectPoint(pnts_init,geom_tol,uvw_list);
        end

        function [uvw_list,unconv_idx]=projectPoint(self,pnts_init,geom_tol,uvw_list)
            % adjust U, V, W by Jacobian transformation
            % also can project point to volume
            %
            % input:
            % pnt_list(matrix): point_number x dimension matrix
            %
            if nargin < 4
                uvw_list=[];
                if nargin < 3
                    geom_tol=[];
                end
            end
            if isempty(geom_tol), geom_tol=sqrt(eps);end
            self=self.deriv(1);

            % find point to start
            if isempty(uvw_list)
                uvw_list=self.findNearest(pnts_init,20);
            end
            [pnt_num,~]=size(pnts_init);
            u_min=min(self.u_knotvctr);u_max=max(self.u_knotvctr);
            v_min=min(self.v_knotvctr);v_max=max(self.v_knotvctr);
            w_min=min(self.w_knotvctr);w_max=max(self.w_knotvctr);
            uvw_min=[u_min,v_min,w_min];uvw_max=[u_max,v_max,w_max];

            % iteration
            iter=0;iter_max=50;
            done=false;
            unconv_idx=(1:pnt_num)';
            while ~done
                [pnts,dpnt_duvs]=self.calGradient(uvw_list(unconv_idx,:));
                dpnts_du=dpnt_duvs{1};
                dpnts_dv=dpnt_duvs{2};
                dpnts_dw=dpnt_duvs{3};
                dpnt_list=pnts_init(unconv_idx,:)-pnts;

                % Jacobian transformation
                RU_RU=sum(dpnts_du.*dpnts_du,2);
                RU_RV=sum(dpnts_du.*dpnts_dv,2);
                RU_RW=sum(dpnts_du.*dpnts_dw,2);
                RV_RV=sum(dpnts_dv.*dpnts_dv,2);
                RV_RW=sum(dpnts_dv.*dpnts_dw,2);
                RW_RW=sum(dpnts_dw.*dpnts_dw,2);
                RU_D=sum(dpnts_du.*dpnt_list,2);
                RV_D=sum(dpnts_dv.*dpnt_list,2);
                RW_D=sum(dpnts_dw.*dpnt_list,2);
                jaco_base=-RW_RW.*RU_RV.^2+2.*RU_RV.*RU_RW.*RV_RW-RV_RV.*RU_RW.^2-RU_RU.*RV_RW.^2+RU_RU.*RV_RV.*RW_RW;
                dus=(-RV_RW.^2.*RU_D+RU_D.*RV_RV.*RW_RW+RV_D.*RU_RW.*RV_RW-RV_D.*RU_RV.*RW_RW+RW_D.*RU_RV.*RV_RW-RW_D.*RU_RW.*RV_RV)./jaco_base;
                dvs=(-RU_RW.^2.*RV_D+RU_D.*RU_RW.*RV_RW-RU_D.*RU_RV.*RW_RW+RV_D.*RU_RU.*RW_RW+RW_D.*RU_RV.*RU_RW-RW_D.*RU_RU.*RV_RW)./jaco_base;
                dws=(-RU_RV.^2.*RW_D+RU_D.*RU_RV.*RV_RW-RU_D.*RU_RW.*RV_RV+RV_D.*RU_RV.*RU_RW-RV_D.*RU_RU.*RV_RW+RW_D.*RU_RU.*RV_RV)./jaco_base;
                dus(isnan(dus) | isinf(dus))=0;
                dvs(isnan(dvs) | isinf(dvs))=0;
                dws(isnan(dws) | isinf(dws))=0;

                uvw_list(unconv_idx,:)=uvw_list(unconv_idx,:)+[dus,dvs,dws];
                uvw_list=max(uvw_list,uvw_min);uvw_list=min(uvw_list,uvw_max);

                unconv_idx=unconv_idx((abs(RU_D)+abs(RV_D)+abs(RW_D) > geom_tol));

                % pnts_inv=self.calPoint(uvw_list);
                % scatter3(pnts_inv(1,:),pnts_inv(2,:),pnts_inv(3,:));

                iter=iter+1;
                if isempty(unconv_idx) || iter >= iter_max
                    done=true;
                end
            end
        end

        function [uvw_list,dis_list]=findNearest(self,pnt_list,param)
            % find nearest U, V, W in grid
            %
            % input:
            % pnt_list(matrix): point_number x dimension matrix
            % param(double): optional, default is 20, grid number to search
            %
            if nargin < 3, param=[];end
            if isempty(param), param=20;end
            pnt_num=size(pnt_list,1);
            uvw_list=zeros(pnt_num,3);
            dis_list=zeros(pnt_num,1);

            % generate rough mesh to initialize pre coord
            base_list=((0:(param-1))+(1:param))/2/param;
            u_base_list=base_list.*(max(self.u_knotvctr)-min(self.u_knotvctr))+min(self.u_knotvctr);
            v_base_list=base_list.*(max(self.v_knotvctr)-min(self.v_knotvctr))+min(self.v_knotvctr);
            w_base_list=base_list.*(max(self.w_knotvctr)-min(self.w_knotvctr))+min(self.w_knotvctr);
            pnt_base_list=self.calPoint({u_base_list,v_base_list,w_base_list});
            [U,V,W]=ndgrid(u_base_list,v_base_list,w_base_list);
            uvw_base_list=[U(:),V(:),W(:)];
            for pnt_idx=1:pnt_num
                pnt=reshape(pnt_list(pnt_idx,:),1,1,1,[]);
                dis_abs=vecnorm(pnt_base_list-pnt,2,4);
                [dis,idx]=min(dis_abs,[],"all");
                uvw_list(pnt_idx,:)=uvw_base_list(idx(1),:);
                dis_list(pnt_idx)=dis;
            end
        end
    end
end
