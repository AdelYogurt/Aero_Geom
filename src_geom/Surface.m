classdef Surface
    % Non-Uniform Rational B-Splines Surface
    % define reference to step standard
    %
    properties % Explicit Attributes
        form='.B-NURBS.'; % surface_form
        coef_dim=[]; % dimension
        u_coef_num; % u number of coefs
        v_coef_num; % v number of coefs
        coefs; % control_points_list
        u_knotvctr; % u_knot_vector
        v_knotvctr; % v_knot_vector
        u_order=[]; % u_degree
        v_order=[]; % v_degree
    end

    properties % Derivate Attributes
        % step standard properties
        knot_spec=''; % knot_spec
        u_closed=false; % UPeriodic (boolean)
        v_closed=false; % VPeriodic (boolean)
        intersected=false; % self_intersect (boolean)

        u_deriv_srf=Surface.empty(0);
        v_deriv_srf=Surface.empty(0);
    end

    methods % define surface
        function self=Surface(poles,u_degree,v_degree,u_mults,v_mults,u_knots,v_knots,weights)
            % generate Non-Uniform Rational B-Splines Surface
            %
            % calling:
            % srf=Surface(poles,u_degree,v_degree,u_mults,v_mults,u_knots,v_knots,weights)
            % srf=Surface(coefs,{u_knots,v_knots})
            %
            % input:
            % poles (matrix): control point, u_pole_num x v_pole_num x dimension matrix
            % u_degree (matrix): optional input
            % v_degree (matrix): optional input
            % u_mults (matrix): optional input
            % v_mults (matrix): optional input
            % u_knots (matrix): optional input
            % v_knots (matrix): optional input
            % weights (matrix): optional input
            %
            % output:
            % Surface
            %
            % notice:
            % degree default is pole_num-1, which will be Bezier surface
            %
            if nargin < 8
                weights=[];
                if nargin < 7
                    v_knots=[];
                    if nargin < 6
                        u_knots=[];
                        if nargin < 5
                            v_mults=[];
                            if nargin < 4
                                u_mults=[];
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

            [u_pole_num,v_pole_num,coef_dim]=size(poles);

            % default value
            if isempty(u_degree),u_degree=u_pole_num-1;end
            if isempty(v_degree),v_degree=v_pole_num-1;end

            % different input mode
            if iscell(u_degree) && length(u_degree) > 1
                coefs=poles;
                [u_knotvrtc,v_knotvrtc]=u_degree{:};
                u_knotvrtc=sort(u_knotvrtc);
                v_knotvrtc=sort(v_knotvrtc);

                u_order=length(u_knotvrtc)-u_pole_num-1;
                v_order=length(v_knotvrtc)-v_pole_num-1;

                if u_pole_num < (u_order+1) || v_pole_num < (v_order+1)
                    error('Surface: pole_num less than order+1');
                end
            elseif length(u_degree) == 1
                % u default value
                if isempty(u_mults) && isempty(u_knots)
                    u_mults=[u_degree+1,ones(1,u_pole_num-u_degree-1),u_degree+1];
                    u_knots=linspace(0,1,u_pole_num-u_degree+1);
                elseif ~isempty(u_mults) && isempty(u_knots)
                    u_knots=linspace(0,1,length(u_mults));
                elseif isempty(u_mults) && ~isempty(u_knots)
                    error('Surface: need u_mults input');
                end
                u_knotvrtc=baseKnotVctr(u_mults,u_knots);
                u_knotvrtc=sort(u_knotvrtc);

                % v default value
                if isempty(v_mults) && isempty(v_knots)
                    v_mults=[v_degree+1,ones(1,v_pole_num-v_degree-1),v_degree+1];
                    v_knots=linspace(0,1,v_pole_num-v_degree+1);
                elseif ~isempty(v_mults) && isempty(v_knots)
                    v_knots=linspace(0,1,length(v_mults));
                elseif isempty(v_mults) && ~isempty(v_knots)
                    error('Surface: need v_mults input');
                end
                v_knotvrtc=baseKnotVctr(v_mults,v_knots);
                v_knotvrtc=sort(v_knotvrtc);

                if isempty(weights), weights=ones(u_pole_num,v_pole_num,1);end
                if ~isempty(weights), weights=reshape(weights,[u_pole_num,v_pole_num,1]);end

                if u_pole_num < (u_degree+1) || v_pole_num < (v_degree+1)
                    error('Surface: pole_num less than degree+1');
                end

                if length(u_knotvrtc) ~= u_pole_num+u_degree+1 || length(v_knotvrtc) ~= v_pole_num+v_degree+1
                    error('Surface: knot_num is not equal to pole_num+degree+1');
                end

                u_order=u_degree;
                v_order=v_degree;
                coefs=cat(3,poles.*weights,weights);
                coef_dim=coef_dim+1;
            else
                error('Surface: error input format');
            end

            self.coef_dim=coef_dim;
            self.u_coef_num=u_pole_num;
            self.v_coef_num=v_pole_num;
            self.u_order=u_order;
            self.v_order=v_order;
            self.coefs=coefs;
            self.u_knotvctr=u_knotvrtc;
            self.v_knotvctr=v_knotvrtc;
        end

        function self=deriv(self,deriv_time)
            % generate derivate BSpline for calculate gradient
            %

            % recursion to derivate data
            if deriv_time > 0
                if isempty(self.u_deriv_srf)
                    % generate derivate surface along the u direction
                    coefs_mat=permute(self.coefs,[1,2,3]);
                    coefs_mat=reshape(coefs_mat,[self.u_coef_num,self.v_coef_num*self.coef_dim]);
                    [coefs_mat,~,u_knotvrtc]=BSpline.deriv(coefs_mat,self.u_order,self.u_knotvctr);
                    coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.v_coef_num,self.coef_dim]);
                    coefs_mat=permute(coefs_mat,[1,2,3]);

                    self.u_deriv_srf=Surface(coefs_mat,{u_knotvrtc,self.v_knotvctr});
                end

                if isempty(self.v_deriv_srf)
                    % generate derivate surface along the v direction
                    coefs_mat=permute(self.coefs,[2,1,3]);
                    coefs_mat=reshape(coefs_mat,[self.v_coef_num,self.u_coef_num*self.coef_dim]);
                    [coefs_mat,~,v_knotvrtc]=BSpline.deriv(coefs_mat,self.v_order,self.v_knotvctr);
                    coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.u_coef_num,self.coef_dim]);
                    coefs_mat=permute(coefs_mat,[2,1,3]);

                    self.v_deriv_srf=Surface(coefs_mat,{self.u_knotvctr,v_knotvrtc});
                end

                % donot output data to avoid matlab repeat deriv_crv
                self.v_deriv_srf=self.v_deriv_srf.deriv(deriv_time-1);
                self.u_deriv_srf=self.u_deriv_srf.deriv(deriv_time-1);
            end
        end

        function [poles,weights]=getPoles(self)
            % convert coefs to poles and weights
            %
            poles=self.coefs(:,:,1:end-1);
            weights=self.coefs(:,:,end);
            weights(weights == 0)=1;
            poles=poles./weights;
        end

        function [pnts,wts]=calPoint(self,us,vs)
            % calculate point on surface
            %
            % calling:
            % [pnts,wts]=srf.calPoint({u_list,v_list})
            % [pnts,wts]=srf.calPoint(uv_list)
            % [pnts,wts]=srf.calPoint(us,vs)
            %
            % notice:
            % pnts and wts is ndgrid generate format
            % last dimension of uvs should equal to 2
            %
            if nargin < 3,vs=[];end
            if isempty(vs),uvs=us;end

            if iscell(us)
                us=uvs{1}(:);vs=uvs{2}(:);
            else
                uvs=cat(3,us,vs);
                siz=size(uvs);uvs=reshape(uvs,[],2);
                us=uvs(:,1);vs=uvs(:,2);
            end

            [u_N_list,u_idx_srt,~]=baseFcnN(us,self.u_order,self.u_knotvctr);
            [v_N_list,v_idx_srt,~]=baseFcnN(vs,self.v_order,self.v_knotvctr);

            if iscell(uvs)
                % input is {u_list,v_list}, mesh point
                u_num=length(us);v_num=length(vs);
                
                % evaluate along the u direction
                coefs_u=reshape(self.coefs,[self.u_coef_num,self.v_coef_num*self.coef_dim]);
                pnts=zeros(length(us),size(coefs_u,2));
                for deg_idx=1:self.u_order+1
                    pnts=pnts+u_N_list(:,deg_idx).*coefs_u(u_idx_srt+(deg_idx-1),:);
                end
                pnts=reshape(pnts,[u_num,self.v_coef_num,self.coef_dim]);

                % evaluate along the v direction
                coefs_v=permute(pnts,[2,1,3]);
                coefs_v=reshape(coefs_v,self.v_coef_num,self.coef_dim*u_num);
                pnts=zeros(length(vs),size(coefs_v,2));
                for deg_idx=1:self.v_order+1
                    pnts=pnts+v_N_list(:,deg_idx).*coefs_v(v_idx_srt+(deg_idx-1),:);
                end
                pnts=reshape(pnts,[v_num,u_num,self.coef_dim]);
                pnts=permute(pnts,[2,1,3]);

                wts=pnts(:,:,end);
                pnts=pnts(:,:,1:end-1);
                if nargout < 2
                    pnts=pnts./wts;
                end
            else
                % input is uv_list, scatter point
                pnt_num=size(uvs,1);

                % evaluate along the u direction
                coefs_u=reshape(self.coefs,[self.u_coef_num,self.v_coef_num*self.coef_dim]);
                pnts=zeros(length(us),size(coefs_u,2));
                for deg_idx=1:self.u_order+1
                    pnts=pnts+u_N_list(:,deg_idx).*coefs_u(u_idx_srt+(deg_idx-1),:);
                end

                % evaluate along the v direction
                coefs_v=reshape(pnts,pnt_num*self.v_coef_num,self.coef_dim);
                pnts=zeros(pnt_num,self.coef_dim);
                for deg_idx=1:self.v_order+1
                    idx=sub2ind([pnt_num,self.v_coef_num],(1:pnt_num)',v_idx_srt+(deg_idx-1));
                    pnts=pnts+v_N_list(:,deg_idx).*coefs_v(idx,:);
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

        function [pnts,us,vs]=calGeom(self,param)
            if nargin < 2,param=[];end
            if isempty(param), param={51,51};end
            if isnumeric(param), param={param,param};end
            u_param=param{1};v_param=param{2};
            if length(u_param) == 1, u_param=linspace(self.u_knotvctr(1),self.u_knotvctr(end),u_param);end
            if length(v_param) == 1, v_param=linspace(self.v_knotvctr(1),self.v_knotvctr(end),v_param);end
            [us,vs]=ndgrid(u_param,v_param);
            pnts=self.calPoint({u_param,v_param});
        end

        function varargout=calGradient(self,us,vs)
            % calculate gradient of surface
            %
            % calling:
            % srf.calPoint({u_list,v_list})
            % srf.calPoint(uv_list)
            % srf.calPoint(us,vs)
            %
            if nargin < 3,vs=[];end
            deriv_time=nargout-1;
            self=self.deriv(deriv_time);

            % recursion to calculate gradient
            if deriv_time > 0
                % calculate initial curve point and weight
                [c0,w0]=self.calPoint(us,vs);
                p0=c0./w0;

                % calculate first derivative
                [c1u,w1u]=self.u_deriv_srf.calPoint(us,vs);
                [c1v,w1v]=self.v_deriv_srf.calPoint(us,vs);
                p1u=(c1u-p0.*w1u)./w0;
                p1v=(c1v-p0.*w1v)./w0;
                p1={p1u,p1v};

                varargout={p0,p1};
                if nargout > 2
                    % calculate second derivative
                    [c2uu,w2uu]=self.u_deriv_srf.u_deriv_srf.calPoint(us,vs);
                    [c2uv,w2uv]=self.u_deriv_srf.v_deriv_srf.calPoint(us,vs);
                    [c2vv,w2vv]=self.v_deriv_srf.v_deriv_srf.calPoint(us,vs);
                    p2uu=((c2uu-2*p1{1}.*w1u)./w0-p0.*w2uu);
                    p2uv=(((c2uv-p1v.*w1u-p1u.*w1v-p0.*w2uv)./w0));
                    p2vv=(c2vv-2*p1{2}.*w1v)./w0-p0.*w2vv;
                    p2={p2uu,p2uv;
                        p2uv,p2vv};
                    varargout=[varargout,{p2}];

                    if nargout > 3
                        % calculate thrid derivative
                        [c3uuu,w3uuu]=self.u_deriv_srf.u_deriv_srf.u_deriv_srf.calPoint(us,vs);
                        [c3uuv,w3uuv]=self.u_deriv_srf.u_deriv_srf.v_deriv_srf.calPoint(us,vs);
                        [c3uvv,w3uvv]=self.u_deriv_srf.v_deriv_srf.v_deriv_srf.calPoint(us,vs);
                        [c3vvv,w3vvv]=self.v_deriv_srf.v_deriv_srf.v_deriv_srf.calPoint(us,vs);
                        p3uuu=(c3uuu-3*p1u.*w2uu-c0.*w3uuu-3*p2uu.*w1u)./w0;
                        p3uuv=(c3uuv-2*p1u.*w2uv-p1v.*w2uu-c0.*w3uuv-2*p2uv.*w1u-p2uu.*w1v)./w0;
                        p3uvu=p3uuv;
                        p3vuu=p3uuv;
                        p3uvv=(c3uvv-2*p1v.*w2uv-p1u.*w2vv-c0.*w3uvv-2*p2uv.*w1v-p2vv.*w1u)./w0;
                        p3vuv=p3uvv;
                        p3vvu=p3uvv;
                        p3vvv=(c3vvv-3*p1v.*w2vv-c0.*w3vvv-3*p2vv.*w1v)./w0;
                        p3=cell(2,2,2);
                        p3{1,1,1}=p3uuu;
                        p3{1,1,2}=p3uuv;
                        p3{1,2,1}=p3uvu;
                        p3{2,1,1}=p3vuu;
                        p3{1,2,2}=p3uvv;
                        p3{2,1,2}=p3vuv;
                        p3{2,2,1}=p3vvu;
                        p3{2,2,2}=p3vvv;
                        varargout=[varargout,{p3}];

                        if nargout > 4
                            % calculate fourth derivative
                            [c4uuuu,w4uuuu]=self.u_deriv_srf.u_deriv_srf.u_deriv_srf.u_deriv_srf.calPoint(us,vs);
                            [c4uuuv,w4uuuv]=self.u_deriv_srf.u_deriv_srf.u_deriv_srf.v_deriv_srf.calPoint(us,vs);
                            [c4uuvv,w4uuvv]=self.u_deriv_srf.u_deriv_srf.v_deriv_srf.v_deriv_srf.calPoint(us,vs);
                            [c4uvvv,w4uvvv]=self.u_deriv_srf.v_deriv_srf.v_deriv_srf.v_deriv_srf.calPoint(us,vs);
                            [p4vvvv,w4vvvv]=self.v_deriv_srf.v_deriv_srf.v_deriv_srf.v_deriv_srf.calPoint(us,vs);
                            p4{1,1,1,1}=(c4uuuu-c0.*w4uuuu-4*p1u.*w3uuu -6*p2uu.*w2uu-4*p3uuu.*w1u)./w0;
                            p4{1,1,1,2}=(c4uuuv-c0.*w4uuuv-3*p1u.*w3uuv-p1v.*w3uuu -3*p2uv.*w2uu -3*p2uu.*w2uv ...
                                -p3uuu.*w1v-3*p3uuv.*w1u)./w0;
                            p4{1,1,2,1}=p4{1,1,1,2};
                            p4{1,2,1,1}=p4{1,1,1,2};
                            p4{2,1,1,1}=p4{1,1,1,2};
                            p4{1,1,2,2}=(c4uuvv-c0.*w4uuvv-2*p1u.*w3uvv-2*p1v.*w3uuv -p2uu.*w2vv -p2vv.*w2uu ...
                                -4*p2uv.*w2uv-2*p3uuv.*w1v-2*p3uvv.*w1u)./w0;
                            p4{1,2,2,1}=p4{1,1,2,2};
                            p4{2,2,1,1}=p4{1,1,2,2};
                            p4{1,2,1,2}=p4{1,1,2,2};
                            p4{2,1,2,1}=p4{1,1,2,2};
                            p4{2,1,1,2}=p4{1,1,2,2};
                            p4{1,2,2,2}=(c4uvvv-c0.*w4uvvv-3*p1v.*w3uvv-p1u.*w3vvv -3*p2uv.*w2vv -3*p2vv.*w2uv ...
                                -p3vvv.*w1u-3*p3uvv.*w1v)./w0;
                            p4{2,2,1,2}=p4{1,2,2,2};
                            p4{2,1,2,2}=p4{1,2,2,2};
                            p4{2,2,2,1}=p4{1,2,2,2};
                            p4{2,2,2,2}=(p4vvvv-c0.*w4vvvv-4*p1v.*w3vvv -6*p2vv.*w2vv-4*p3vvv.*w1v)./w0;
                            varargout=[varargout,{p4}];
                        end
                    end
                end
            end
        end

        function varargout=displayGeom(self,axe_hdl,srf_option,param)
            % draw surface on axes handle
            %
            if nargin < 4
                param=[];
                if nargin < 3
                    srf_option=[];
                    if nargin < 2
                        axe_hdl=[];
                    end
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end
            if isempty(srf_option), srf_option=struct('LineStyle','none','FaceAlpha',0.5);end
            if isempty(param), param={51,51};end
            if isnumeric(param), param={param,param};end
            u_param=param{1};v_param=param{2};
            if length(u_param) == 1, u_param=linspace(self.u_knotvctr(1),self.u_knotvctr(end),u_param);end
            if length(v_param) == 1, v_param=linspace(self.v_knotvctr(1),self.v_knotvctr(end),v_param);end
            pnts=self.calPoint({u_param,v_param});

            % draw points on axe_hdl
            if self.coef_dim-1 == 2
                srf_hdl=surface(axe_hdl,pnts(:,:,1)',pnts(:,:,2)',zeros(size(pnts(:,:,1)')),srf_option);
            else
                srf_hdl=surface(axe_hdl,pnts(:,:,1)',pnts(:,:,2)',pnts(:,:,3)',srf_option);
                zlabel('\itZ');
            end
            xlabel('\itX');
            ylabel('\itY');

            if nargout > 0, varargout={srf_hdl};end
        end

        function varargout=displayPole(self,axe_hdl,pole_option)
            % draw curve on figure handle
            %
            if nargin < 3
                pole_option=[];
                if nargin < 2
                    axe_hdl=[];
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % default draw option
            if isempty(pole_option)
                pole_option=struct('Marker','s','MarkerEdgeColor','r','EdgeColor','r','LineStyle','--','FaceAlpha',0);
            end

            % draw poles on axe_hdl
            poles=self.getPoles();
            if self.coef_dim-1 == 2
                srf_hdl=surface(axe_hdl,poles(:,:,1)',poles(:,:,2)',pole_option);
            else
                srf_hdl=surface(axe_hdl,poles(:,:,1)',poles(:,:,2)',poles(:,:,3)',pole_option);
                zlabel('\itZ');
            end
            xlabel('\itX');
            ylabel('\itY');

            if nargout > 0, varargout={srf_hdl};end
        end

        function varargout=displayDirect(self,axe_hdl)
            % draw curve on figure handle
            %
            if nargin < 2
                axe_hdl=[];
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % calculate coord on curve
            pnt_mesh=self.coefs(1:2,1:2,1:end-1);
            wt_mesh=self.coefs(1:2,1:2,end);
            wt_mesh(wt_mesh == 0)=1;
            pnt_mesh=pnt_mesh./wt_mesh;
            origin=pnt_mesh(1,1,:);
            coord=squeeze([pnt_mesh(2,1,:)-origin,pnt_mesh(1,2,:)-origin]);
            coord=(coord./vecnorm(coord,2,2))';

            % scale of quiver
            scale=norm(diff([axe_hdl.XLim;axe_hdl.YLim;axe_hdl.ZLim],1,2));
            coord=coord*scale*0.05;

            % draw coord on axe_hdl
            if self.coef_dim-1 == 2
                hold on;
                quv_hdl(1)=quiver(axe_hdl,origin(1),origin(2),coord(1,1),coord(2,1),0,'Color','r');
                quv_hdl(2)=quiver(axe_hdl,origin(1),origin(2),coord(1,2),coord(2,2),0,'Color','g');
                hold off;
            elseif self.coef_dim-1 == 3
                hold on;
                quv_hdl(1)=quiver3(axe_hdl,origin(1),origin(2),origin(3),coord(1,1),coord(2,1),coord(3,1),0,'Color','r');
                quv_hdl(2)=quiver3(axe_hdl,origin(1),origin(2),origin(3),coord(1,2),coord(2,2),coord(3,2),0,'Color','g');
                hold off;
                zlabel('\itZ');
                
            end
            xlabel('\itX');
            ylabel('\itY');

            if nargout > 0, varargout={quv_hdl};end
        end

        function varargout=extract(self,sides)
            % get boundary curve
            %
            if nargin < 2,sides=[];end
            if isempty(sides), sides=1:4;end

            crv_list=Curve.empty(0,4);
            for idx=1:2
                srf_idx1=(idx-1)*2+1;
                srf_idx2=(idx-1)*2+2;

                switch idx
                    case 1 % U
                        coefs1=squeeze(self.coefs(1,:,:));
                        coefs2=squeeze(self.coefs(end,:,:));
                        crv_list(srf_idx1)=Curve(coefs1,self.v_knotvctr);
                        crv_list(srf_idx2)=Curve(coefs2,self.v_knotvctr);
                    case 2 % V
                        coefs1=squeeze(self.coefs(:,1,:));
                        coefs2=squeeze(self.coefs(:,end,:));
                        crv_list(srf_idx1)=Curve(coefs1,self.u_knotvctr);
                        crv_list(srf_idx2)=Curve(coefs2,self.u_knotvctr);
                end
            end

            crv_list=crv_list(sides);
            if nargout > 1, varargout=num2cell(crv_list);
            else, varargout={crv_list};end
        end
    end

    methods % control surface
        function self=permute(self,ord)
            % change surface parameter order
            %
            knots_list={self.u_knotvctr,self.v_knotvctr};
            self=Surface(permute(self.coefs,[ord,3]),knots_list(ord));
        end

        function self=reverseU(self)
            % revese U direction of surface
            %
            self.coefs=self.coefs(end:-1:1,:,:);
            min_k=min(self.u_knotvctr);max_k=max(self.u_knotvctr);dk=max_k-min_k;
            self.u_knotvctr=dk-(self.u_knotvctr(end:-1:1)-min_k)+min_k;
        end

        function self=reverseV(self)
            % revese V direction of surface
            %
            self.coefs=self.coefs(:,end:-1:1,:);
            min_k=min(self.v_knotvctr);max_k=max(self.v_knotvctr);dk=max_k-min_k;
            self.v_knotvctr=dk-(self.v_knotvctr(end:-1:1)-min_k)+min_k;
        end

        function self=addOrder(self,u_times,v_times)
            % increase surface degree
            %
            if u_times <= 0 && v_times <= 0, return;end

            % modify order along the u direction
            if u_times <= 0
                coefs_mat=permute(self.coefs,[1,2,3]);
                coefs_mat=reshape(coefs_mat,[self.u_coef_num,self.v_coef_num*self.coef_dim]);
                [coefs_mat,self.u_order,self.u_knotvctr]=BSpline.addOrder(coefs_mat,self.u_order,self.u_knotvctr,u_times);
                coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.v_coef_num,self.coef_dim]);
                self.coefs=permute(coefs_mat,[1,2,3]);
                self.u_coef_num=size(self.coefs,1);
            end

            % modify order along the v direction
            if v_times <= 0
                coefs_mat=permute(self.coefs,[2,1,3]);
                coefs_mat=reshape(coefs_mat,[self.v_coef_num,self.u_coef_num*self.coef_dim]);
                [coefs_mat,self.v_order,self.v_knotvctr]=BSpline.addOrder(coefs_mat,self.v_order,self.v_knotvctr,v_times);
                coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.u_coef_num,self.coef_dim]);
                self.coefs=permute(coefs_mat,[2,1,3]);
                self.v_coef_num=size(self.coefs,2);
            end
        end

        function self=insertKnot(self,u_knot_ins,v_knot_ins)
            % insert knot to surface
            %
            if isempty(u_knot_ins) && isempty(v_knot_ins), return;end

            % insert knots along the u direction
            if ~isempty(u_knot_ins)
                coefs_mat=permute(self.coefs,[1,2,3]);
                coefs_mat=reshape(coefs_mat,[self.u_coef_num,self.v_coef_num*self.coef_dim]);
                [coefs_mat,self.u_order,self.u_knotvctr]=BSpline.insertKnot(coefs_mat,self.u_order,self.u_knotvctr,u_knot_ins);
                coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.v_coef_num,self.coef_dim]);
                self.coefs=permute(coefs_mat,[1,2,3]);
                self.u_coef_num=size(self.coefs,1);
            end

            % insert knots along the v direction
            if ~isempty(v_knot_ins)
                coefs_mat=permute(self.coefs,[2,1,3]);
                coefs_mat=reshape(coefs_mat,[self.v_coef_num,self.u_coef_num*self.coef_dim]);
                [coefs_mat,self.v_order,self.v_knotvctr]=BSpline.insertKnot(coefs_mat,self.v_order,self.v_knotvctr,v_knot_ins);
                coefs_mat=reshape(coefs_mat,[size(coefs_mat,1),self.u_coef_num,self.coef_dim]);
                self.coefs=permute(coefs_mat,[2,1,3]);
                self.v_coef_num=size(self.coefs,2);
            end
        end

        function self=translate(self,tran_vctr)
            % translate surface
            %
            tran_vctr=reshape(tran_vctr,1,1,[]);
            self.coefs(:,:,1:end-1)=self.coefs(:,:,1:end-1)+self.coefs(:,:,end).*tran_vctr(1:self.coef_dim-1);
        end

        function self=rotate(self,rot_mat,rot_cntr)
            % rotate surface
            %
            if nargin < 3 || isempty(rot_cntr), rot_cntr=zeros(self.coef_dim-1,1);end
            self=self.translate(-rot_cntr);
            self.coefs(:,:,1:end-1)=reshape(reshape(self.coefs(:,:,1:end-1),[],self.coef_dim-1)*rot_mat',self.u_coef_num,self.v_coef_num,self.coef_dim-1);
            self=self.translate(rot_cntr);
        end
    end

    methods % calculate coord
        function uv_list=calCoordinate(self,pnts_init,geom_tol)
            % base on X, Y, Z calculate local coordinate in surface
            %
            if nargin < 3, geom_tol=[];end
            if isempty(geom_tol), geom_tol=sqrt(eps);end

            % find point to start
            uv_list=self.findNearest(pnts_init,20);

            % use project function to adjust parameter
            uv_list=self.projectPoint(pnts_init,geom_tol,uv_list);
        end

        function uv_list=projectPoint(self,pnts_init,geom_tol,uv_list)
            % adjust U, V by Jacobian transformation
            % also can project point to surface
            %
            % input:
            % pnt_list(matrix): point_number x dimension matrix
            %
            if nargin < 4
                uv_list=[];
                if nargin < 3
                    geom_tol=[];
                end
            end
            if isempty(geom_tol), geom_tol=sqrt(eps);end
            self=self.deriv(1);

            % find point to start
            if isempty(uv_list)
                uv_list=self.findNearest(pnts_init,20);
            end
            [pnt_num,~]=size(pnts_init);
            uv_min=[self.u_knotvctr(1),self.v_knotvctr(1)];
            uv_max=[self.u_knotvctr(end),self.v_knotvctr(end)];

            % iteration
            iter=0;iter_max=50;
            done=false;
            pnt_idx=(1:pnt_num)';
            while ~done
                [pnts,dpnt_duvs]=self.calGradient(uv_list(pnt_idx,:));
                dpnts_du=dpnt_duvs{1};
                dpnts_dv=dpnt_duvs{2};
                dpnts=pnts_init(pnt_idx,:)-pnts;

                % Jacobian transformation
                RU_RU=sum(dpnts_du.*dpnts_du,2);
                RU_RV=sum(dpnts_du.*dpnts_dv,2);
                RV_RV=sum(dpnts_dv.*dpnts_dv,2);
                RU_D=sum(dpnts_du.*dpnts,2);
                RV_D=sum(dpnts_dv.*dpnts,2);
                jaco_base=RU_RU.*RV_RV-RU_RV.*RU_RV;
                dus=(RU_D.*RV_RV-RV_D.*RU_RV)./jaco_base;
                dvs=(RV_D.*RU_RU-RU_D.*RU_RV)./jaco_base;
                dus(isnan(dus) | isinf(dus))=0;
                dvs(isnan(dvs) | isinf(dvs))=0;

                uv_list(pnt_idx,:)=uv_list(pnt_idx,:)+[dus,dvs];
                uv_list=max(uv_list,uv_min);uv_list=min(uv_list,uv_max);

                pnt_idx=pnt_idx((abs(RU_D)+abs(RV_D) > geom_tol));

                % pnts_inv=self.calPoint(uv_list);
                % scatter3(pnts_inv(:,1),pnts_inv(:,2),pnts_inv(:,3));

                iter=iter+1;
                if isempty(pnt_idx) || iter >= iter_max
                    done=true;
                end
            end
        end

        function uv_list=findNearest(self,pnt_list,param)
            % find nearest U, V in grid
            %
            % input:
            % pnt_list(matrix): point_number x dimension matrix
            % param(double): optional, default is 20, grid number to search
            %
            if nargin < 3, param=[];end
            if isempty(param), param=20;end
            pnt_num=size(pnt_list,1);

            % generate rough mesh to initialize pre coord
            base_list=((0:(param-1))+(1:param))/2/param;
            u_base_list=base_list*(self.u_knotvctr(end)-self.u_knotvctr(1))+self.u_knotvctr(1);
            v_base_list=base_list*(self.v_knotvctr(end)-self.v_knotvctr(1))+self.v_knotvctr(1);
            pnt_base_list=self.calPoint({u_base_list,v_base_list});
            [U,V]=ndgrid(u_base_list,v_base_list);
            uv_base_list=[U(:),V(:)];
            uv_list=zeros(pnt_num,2);
            for pnt_idx=1:pnt_num
                pnt=reshape(pnt_list(pnt_idx,:),1,1,[]);
                dis_abs=vecnorm(pnt_base_list-pnt,2,3);
                [~,idx]=min(dis_abs,[],"all");
                uv_list(pnt_idx,:)=uv_base_list(idx(1),:);
            end
        end
    end
end
