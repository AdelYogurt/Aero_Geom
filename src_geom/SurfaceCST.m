classdef SurfaceCST
    % Class-Shape Transformation (CST) Surface
    %
    % notice:
    % origin reference of shape function is local coordinate
    % point calculate by class and shape will be translate to global
    %
    % reference:
    % [1] Kulfan B. A Universal Parametric Geometry Representation Method -
    % "CST" [C]. 45th AIAA Aerospace Sciences Meeting and Exhibit, Reno,
    % Nevada, U.S.A.
    % [2] Liu C, Duan Y, Cai J, et al. Applications of multi-block cst
    % method for quasi-waverider design[C]. 52nd Aerospace Sciences
    % Meeting, National Harbor, Maryland, U.S.A.
    %
    properties % CST parameter
        form='.CST.'; % curve_form

        % class parameter
        N1x=[];
        N2x=[];
        N1y=[];
        N2y=[];
        N1zv=[];
        N2zv=[];
        N1zu=[];
        N2zu=[];

        % shape parameter
        coef_dim=4; % surface dimension
        u_coef_num=2; % surface u number of coefs
        v_coef_num=2; % surface v number of coefs
        coefs=cat(3,[0,0;1,1],[0,1;0,1],[1,1;1,1]); % surface control_points_list
        u_knotvctr=[0,0,1,1]; % surface u_knot_vector
        v_knotvctr=[0,0,1,1]; % surface v_knot_vector
        u_order=1; % surface u_degree
        v_order=1; % surface v_degree

        shape_fcn=[]; % function handle, pnts=shape_fcn(us, vs)

        % bias parameter
        B11=0;
        B21=0;
        B12=0;
        B22=0;

        bias_fcn=[]; % function handle, pnts=shape_fcn(us, vs)

        coord=[]; % coordination of curve, [Ax1,Ax2,Ax3]
        origin=[]; % origin point of coordination
        sym_x; % if true, u_class fcn will start from 0.5 to 1
        sym_y; % if true, v_class fcn will start from 0.5 to 1
        u_closed=false; % UPeriodic (boolean)
        v_closed=false; % VPeriodic (boolean)
        intersected=false; % self_intersect (boolean)

        u_deriv_srf=[];
        v_deriv_srf=[];
    end

    properties(Access=private) % fit point set parameter
        nodes; % fit point set
        u_nodes; % u fit point parameter
        v_nodes; % v fit point parameter
        u_fit_mat; % u fit matrix
        v_fit_mat; % v fit matrix
    end

    methods % define surface
        function self=SurfaceCST(C_par_x,C_par_y,C_par_zv,C_par_zu,S_par,B_par,sym_x,sym_y)
            % generate CST surface
            %
            % input:
            % C_par_x (list): N1x and N2x array([N1x, N2x]), default is [0.0,0.0]
            % C_par_y (list): N1y and N2y array([N1y, N2y]), default is [0.0,0.0]
            % C_par_zv (list): N1zv and N2zv array([N1zv, N2zv]), default is [0.0,0.0]
            % C_par_zu (list): N1zu and N2zu array([N1zu, N2zu]), default is [0.0,0.0]
            % S_par (list): LX, LY and LZ array([LX, LY, LZ]), default is [1.0,1.0,1.0]
            % B_par (matrix): B11, B21, B12 and B22 matrix([B11, B21; B12, B22]), default is [0.0,0.0;0.0,0.0]
            % sym_x (boolean): whether using half of u while calculate class_fcn
            % sym_y (boolean): whether using half of v while calculate class_fcn
            %
            % output:
            % SurfaceCST
            %
            % notice:
            % if C_par is [0.0, 0.0] or empty, class_fcn will equal to 1.0
            %
            if nargin < 8
                sym_y=[];
                if nargin < 7
                    sym_x=[];
                    if nargin < 6
                        B_par=[];
                        if nargin < 5
                            S_par=[];
                            if nargin < 4
                                C_par_zu=[];
                                if nargin < 3
                                    C_par_zv=[];
                                    if nargin < 2
                                        C_par_y=[];
                                        if nargin < 1
                                            C_par_x=[];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            % default value
            if isempty(sym_y), sym_y=false;end
            if isempty(sym_x), sym_x=false;end
            if isempty(B_par),B_par=[0.0,0.0;0.0,0.0];end
            if isempty(S_par),S_par=[1.0,1.0,1.0];end
            if isempty(C_par_zu),C_par_zu=[0,0];end
            if isempty(C_par_zv),C_par_zv=[0,0];end
            if isempty(C_par_y),C_par_y=[0,0];end
            if isempty(C_par_x),C_par_x=[0,0];end

            self.sym_x=sym_x;
            self.sym_y=sym_y;

            % process C_par

            % process C_par
            self.N1x=C_par_x(1);
            self.N2x=C_par_x(2);
            self.N1y=C_par_y(1);
            self.N2y=C_par_y(2);
            self.N1zv=C_par_zv(1);
            self.N2zv=C_par_zv(2);
            self.N1zu=C_par_zu(1);
            self.N2zu=C_par_zu(2);

            % process S_par
            if isa(S_par,'function_handle')
                shape_fcn=S_par;
                self.shape_fcn=shape_fcn;
                pnt=shape_fcn(0,0);
                coef_dim=length(pnt)+1;

                self.coef_dim=coef_dim;
            else
                LX=S_par(1);LY=S_par(2);LZ=S_par(3);

                pnt_x=[0,0;LX,LX];
                pnt_y=[0,LY;0,LY];
                pnt_z=[LZ,LZ;LZ,LZ];
                pnt_w=ones(2);
                coef_dim=4;u_coef_num=2;v_coef_num=2;u_order=1;v_order=1;
                coefs=cat(3,pnt_x,pnt_y,pnt_z,pnt_w);

                self.coef_dim=coef_dim;
                self.u_coef_num=u_coef_num;
                self.v_coef_num=v_coef_num;
                self.coefs=coefs;
                self.u_order=u_order;
                self.v_order=v_order;
            end

            % process B_par
            if isa(B_par,'function_handle')
                self.bias_fcn=B_par;
            else
                B_par=ones(2).*B_par;
                self.B11=B_par(1,1);
                self.B21=B_par(2,1);
                self.B12=B_par(1,2);
                self.B22=B_par(2,2);
            end

            self.coord=eye(coef_dim-1);
            self.origin=zeros(1,coef_dim-1);
        end

        function self=deriv(self,~)
            % generate derivation data
            %
            if isempty(self.shape_fcn)
                if isempty(self.u_deriv_srf)
                    % generate derivate surface along the u direction
                    dcoefs_du=permute(self.coefs,[1,2,3]);
                    dcoefs_du=reshape(dcoefs_du,[self.u_coef_num,self.v_coef_num*self.coef_dim]);
                    if self.u_coef_num == 2
                        dcoefs_du=dcoefs_du(2,:)-dcoefs_du(1,:);
                        dorder_du=0;dknots_du=[0,1];
                    else
                        [dcoefs_du,dorder_du,dknots_du]=BSpline.deriv(dcoefs_du,self.u_order,self.u_knotvctr);
                    end
                    dcoefs_du=reshape(dcoefs_du,[size(dcoefs_du,1),self.v_coef_num,self.coef_dim]);
                    dcoefs_du=permute(dcoefs_du,[1,2,3]);

                    self.u_deriv_srf.coef_dim=self.coef_dim;
                    self.u_deriv_srf.u_coef_num=self.u_coef_num-1;
                    self.u_deriv_srf.v_coef_num=self.v_coef_num;
                    self.u_deriv_srf.coefs=dcoefs_du;
                    self.u_deriv_srf.u_knots=dknots_du;
                    self.u_deriv_srf.v_knots=self.v_knotvctr;
                    self.u_deriv_srf.u_order=dorder_du;
                    self.u_deriv_srf.v_order=self.v_order;
                end

                if isempty(self.v_deriv_srf)
                    dcoefs_dv=permute(self.coefs,[2,1,3]);
                    dcoefs_dv=reshape(dcoefs_dv,[self.v_coef_num,self.u_coef_num*self.coef_dim]);
                    if self.v_coef_num == 2
                        dcoefs_dv=dcoefs_dv(2,:)-dcoefs_dv(1,:);
                        dorder_dv=0;dknots_dv=[0,1];
                    else
                        [dcoefs_dv,dorder_dv,dknots_dv]=BSpline.deriv(dcoefs_dv,self.v_order,self.v_knotvctr);
                    end
                    dcoefs_dv=reshape(dcoefs_dv,[size(dcoefs_dv,1),self.u_coef_num,self.coef_dim]);
                    dcoefs_dv=permute(dcoefs_dv,[2,1,3]);

                    self.v_deriv_srf.coef_dim=self.coef_dim;
                    self.v_deriv_srf.u_coef_num=self.u_coef_num;
                    self.v_deriv_srf.v_coef_num=self.v_coef_num-1;
                    self.v_deriv_srf.coefs=dcoefs_dv;
                    self.v_deriv_srf.u_knots=self.u_knotvctr;
                    self.v_deriv_srf.v_knots=dknots_dv;
                    self.v_deriv_srf.u_order=self.u_order;
                    self.v_deriv_srf.v_order=dorder_dv;
                end
            end
        end

        function self=addSpline(self,poles,u_degree,v_degree,u_mults,v_mults,u_knots,v_knots,weights)
            % add Bezier/B-Spline surface as shape function
            %
            % calling:
            % srf=srf.addSpline(poles,u_degree,v_degree,u_mults,v_mults,u_knots,v_knots,weights)
            % srf=srf.addSpline(coefs,{u_knots,v_knots})
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
            % SurfaceCST
            %
            % notice:
            % degree default is pole_num-1, which will be Bezier surface
            %
            if nargin < 9
                weights=[];
                if nargin < 8
                    v_knots=[];
                    if nargin < 7
                        u_knots=[];
                        if nargin < 6
                            v_mults=[];
                            if nargin < 5
                                u_mults=[];
                                if nargin < 4
                                    v_degree=[];
                                    if nargin < 3
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
            else
                % u default value
                if isempty(u_mults) && isempty(u_knots)
                    u_mults=[u_degree+1,ones(1,u_pole_num-u_degree-1),u_degree+1];
                    u_knots=linspace(0,1,u_pole_num-u_degree+1);
                elseif ~isempty(u_mults) && isempty(u_knots)
                    u_knots=linspace(0,1,length(u_mults));
                elseif isempty(u_mults) && ~isempty(u_knots)
                    error('SurfaceCST.addSpline: need u_mults input');
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
                    error('SurfaceCST.addSpline: need v_mults input');
                end
                v_knotvrtc=baseKnotVctr(v_mults,v_knots);
                v_knotvrtc=sort(v_knotvrtc);

                if isempty(weights), weights=ones(u_pole_num,v_pole_num,1);end
                if ~isempty(weights), weights=reshape(weights,[u_pole_num,v_pole_num,1]);end

                if u_pole_num < (u_degree+1) || v_pole_num < (v_degree+1)
                    error('SurfaceCST.addSpline: pole_num less than degree+1');
                end

                if length(u_knotvrtc) ~= u_pole_num+u_degree+1 || length(v_knotvrtc) ~= v_pole_num+v_degree+1
                    error('SurfaceCST.addSpline: knot_num is not equal to pole_num+degree+1');
                end

                u_order=u_degree;
                v_order=v_degree;
                coefs=cat(3,poles.*weights,weights);
                coef_dim=coef_dim+1;
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

        function [poles,weights]=getPoles(self)
            % convert coefs to poles and weights
            %
            poles=self.coefs(:,:,1:end-1);
            weights=self.coefs(:,:,end);
            weights(weights == 0)=1;
            poles=poles./weights;
        end

        function [pnts,cls,sps,bias]=calPoint(self,us,vs)
            % calculate point on surface
            %
            % calling:
            % [pnts,wts]=srf.calPoint({u_list,v_list})
            % [pnts,wts]=srf.calPoint(uv_list)
            % [pnts,wts]=srf.calPoint(us,vs)
            %
            % notice:
            % pnts cls, sps and bias is ndgrid generate format, [u,v]
            % last dimension of uvs should equal to 2
            %
            if nargin < 3,vs=[];end
            
            % calculate shape
            if isempty(self.shape_fcn)
                sps=calSurface(self.coef_dim,self.u_coef_num,self.v_coef_num,self.coefs,self.u_knotvctr,self.v_knotvctr,self.u_order,self.v_order,us,vs);
                sps=reshape(sps,[],self.coef_dim-1);
            end

            % calculate perpare
            if isempty(vs),uvs=us;end
            if iscell(us)
                us=uvs{1}(:);vs=uvs{2}(:);
                siz=[length(us),length(vs),self.coef_dim];
                [us,vs]=ndgrid(us,vs);
                us=us(:);vs=vs(:);
            else
                uvs=cat(3,us,vs);
                siz=size(uvs);uvs=reshape(uvs,[],2);
                us=uvs(:,1);vs=uvs(:,2);
            end

            % calculate shape
            if ~isempty(self.shape_fcn)
                sps=reshape(self.shape_fcn(us,vs),[],self.coef_dim-1);
            end

            u_class=us;if self.sym_x,u_class=(u_class/2)+0.5;end
            v_class=vs;if self.sym_y,v_class=(v_class/2)+0.5;end

            % calculate class
            cls_x=baseFcnClass(v_class,self.N1x,self.N2x);
            cls_y=baseFcnClass(u_class,self.N1y,self.N2y);
            cls_zv=baseFcnClass(v_class,self.N1zv,self.N2zv);
            cls_zu=baseFcnClass(u_class,self.N1zu,self.N2zu);
            cls=[cls_x,cls_y,cls_zv.*cls_zu,ones(size(cls_x,1),self.coef_dim-4)];

            % calculate bias
            if isempty(self.bias_fcn)
                bias=self.B11*(1-us).*(1-vs)+self.B21*us.*(1-vs)+self.B12*(1-us).*vs+self.B22*us.*vs;
                bias=[zeros(size(bias,1),2),bias,zeros(size(bias,1),self.coef_dim-4)];
            else
                bias=reshape(self.bias_fcn(us,vs),[],1);
                bias=[zeros(size(bias,1),2),bias,zeros(size(bias,1),self.coef_dim-4)];
            end

            % calculate point
            pnts=cls.*sps+bias;
            pnts=self.convertLocalToGlobal(pnts);

            pnts=reshape(pnts,[siz(1:end-1),self.coef_dim-1]);
            cls=reshape(cls,[siz(1:end-1),self.coef_dim-1]);
            sps=reshape(sps,[siz(1:end-1),self.coef_dim-1]);
            bias=reshape(bias,[siz(1:end-1),self.coef_dim-1]);
        end

        function [pnts,dpnts_duv]=calGradient(self,us,vs,step)
            % calculate gradient of surface
            %
            % calling:
            % srf.calPoint({u_list,v_list})
            % srf.calPoint(uv_list)
            % srf.calPoint(us,vs)
            %
            if nargin < 4
                step=[];
                if nargin < 3
                    vs=[];
                end
            end
            if isempty(step),step=sqrt(eps);end
            self=self.deriv(1);

            % calculate shape
            if isempty(self.shape_fcn)
                [sps_p,sps_w]=calSurface(self.coef_dim,self.u_coef_num,self.v_coef_num,self.coefs,self.u_knotvctr,self.v_knotvctr,self.u_order,self.v_order,us,vs);
                dsrf_du=self.u_deriv_srf;
                dsrf_dv=self.v_deriv_srf;
                [dsps_du_p,dsps_du_w]=calSurface(dsrf_du.coef_dim,dsrf_du.u_coef_num,dsrf_du.v_coef_num,dsrf_du.coefs,dsrf_du.u_knots,dsrf_du.v_knots,dsrf_du.u_order,dsrf_du.v_order,us,vs);
                [dsps_dv_p,dsps_dv_w]=calSurface(dsrf_dv.coef_dim,dsrf_dv.u_coef_num,dsrf_dv.v_coef_num,dsrf_dv.coefs,dsrf_dv.u_knots,dsrf_dv.v_knots,dsrf_dv.u_order,dsrf_dv.v_order,us,vs);
                
                sps=sps_p./sps_w;
                dsps_du=(dsps_du_p-sps.*dsps_du_w)./sps_w;
                dsps_dv=(dsps_dv_p-sps.*dsps_dv_w)./sps_w;

                sps=reshape(sps,[],self.coef_dim-1);
                dsps_du=reshape(dsps_du,[],self.coef_dim-1);
                dsps_dv=reshape(dsps_dv,[],self.coef_dim-1);

                du=self.u_knotvctr(end)-self.u_knotvctr(1);
                if self.u_coef_num ~= 2,dsps_du=dsps_du*du;end
                dv=self.v_knotvctr(end)-self.v_knotvctr(1);
                if self.v_coef_num ~= 2,dsps_dv=dsps_dv*dv;end
            end

            % calculate perpare
            if isempty(vs),uvs=us;end
            if iscell(us)
                us=uvs{1}(:);vs=uvs{2}(:);
                siz=[length(us),length(vs),self.coef_dim];
                [us,vs]=ndgrid(us,vs);
                us=us(:);vs=vs(:);
            else
                uvs=cat(3,us,vs);
                siz=size(uvs);uvs=reshape(uvs,[],2);
                us=uvs(:,1);vs=uvs(:,2);
            end

            % calculate shapes
            if ~isempty(self.shape_fcn)  
                [sps,dsps_du,dsps_dv]=differFcn2Robust(self.shape_fcn,us,vs,step,self.coef_dim-1);
            end

            u_class=us;if self.sym_x,u_class=(u_class/2)+0.5;end
            v_class=vs;if self.sym_y,v_class=(v_class/2)+0.5;end

            % calculate class
            [cls_x,dcls_x_dv]=baseFcnClass(v_class,self.N1x,self.N2x);
            [cls_y,dcls_y_du]=baseFcnClass(u_class,self.N1y,self.N2y);
            [cls_zv,dcls_zv_dv]=baseFcnClass(v_class,self.N1zv,self.N2zv);
            [cls_zu,dcls_zu_du]=baseFcnClass(u_class,self.N1zu,self.N2zu);
            if self.sym_y,dcls_x_dv=dcls_x_dv/2;end,dcls_x_dv=min(dcls_x_dv,1e6);
            if self.sym_x,dcls_y_du=dcls_y_du/2;end,dcls_y_du=min(dcls_y_du,1e6);
            if self.sym_y,dcls_zv_dv=dcls_zv_dv/2;end,dcls_zv_dv=min(dcls_zv_dv,1e6);
            if self.sym_x,dcls_zu_du=dcls_zu_du/2;end,dcls_zu_du=min(dcls_zu_du,1e6);
            cls=[cls_x,cls_y,cls_zv.*cls_zu,ones(size(cls_x,1),self.coef_dim-4)];
            dcls_du=[zeros(size(cls_x,1),1),dcls_y_du,cls_zv.*dcls_zu_du,zeros(size(cls_x,1),self.coef_dim-4)];
            dcls_dv=[dcls_x_dv,zeros(size(cls_y,1),1),dcls_zv_dv.*cls_zu,zeros(size(cls_x,1),self.coef_dim-4)];

            % calculate bias
            if isempty(self.bias_fcn)
                bias=self.B11*(1-us).*(1-vs)+self.B21*us.*(1-vs)+self.B12*(1-us).*vs+self.B22*us.*vs;
                dbias_du=(self.B21-self.B11).*(1-vs)+(self.B22-self.B12).*vs;
                dbias_dv=(self.B12-self.B11).*(1-us)+(self.B22-self.B21).*us;
            else
                [bias,dbias_du,dbias_dv]=differFcn2Robust(self.bias_fcn,us,vs,step,1);
            end
            bias=[zeros(size(bias,1),2),bias,zeros(size(bias,1),self.coef_dim-4)];
            dbias_du=[zeros(size(dbias_du,1),2),dbias_du,zeros(size(dbias_du,1),self.coef_dim-4)];
            dbias_dv=[zeros(size(dbias_dv,1),2),dbias_dv,zeros(size(dbias_dv,1),self.coef_dim-4)];

            pnts=cls.*sps+bias;
            dpnts_du=dcls_du.*sps+cls.*dsps_du+dbias_du;
            dpnts_dv=dcls_dv.*sps+cls.*dsps_dv+dbias_dv;
            pnts=self.convertLocalToGlobal(pnts);
            dpnts_du=dpnts_du*self.coord';
            dpnts_dv=dpnts_dv*self.coord';
            
            pnts=reshape(pnts,[siz(1:end-1),self.coef_dim-1]);
            dpnts_du=reshape(dpnts_du,[siz(1:end-1),self.coef_dim-1]);
            dpnts_dv=reshape(dpnts_dv,[siz(1:end-1),self.coef_dim-1]);

            dpnts_duv={dpnts_du,dpnts_dv};
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
            if length(u_param) == 1, u_param=linspace(0,1,u_param);end
            if length(v_param) == 1, v_param=linspace(0,1,v_param);end
            msh_1zu=max(min(1/self.N1zu,2),0.5);msh_2zu=max(min(1/self.N2zu,2),0.5);
            if (self.N1zu ~= 0) && (self.N2zu ~= 0) && (~self.sym_x),u_param=u_param.^msh_1zu.*(1-u_param)+u_param.*(1-(1-u_param).^msh_2zu);
            elseif (self.N1zu == 0) && (self.N2zu ~= 0) || (self.sym_x),u_param=1-(1-u_param).^(msh_2zu);
            elseif (self.N1zu ~= 0) && (self.N2zu == 0),u_param=u_param.^(msh_1zu);
            end
            msh_1zv=max(min(1/self.N1zv,2),0.5);msh_2zv=max(min(1/self.N2zv,2),0.5);
            if (self.N1zv ~= 0) && (self.N2zv ~= 0) && (~self.sym_y),v_param=v_param.^msh_1zv.*(1-v_param)+v_param.*(1-(1-v_param).^msh_2zv);
            elseif (self.N1zv == 0) && (self.N2zv ~= 0) || (self.sym_y),v_param=1-(1-v_param).^(msh_2zv);
            elseif (self.N1zv ~= 0) && (self.N2zv == 0),v_param=v_param.^(msh_1zv);
            end
            pnts=self.calPoint({u_param,v_param});

            % draw points on axe_hdl
            % if self.coef_dim-1 == 2
            %     srf_hdl=surface(axe_hdl,pnts(:,:,1)',pnts(:,:,2)',srf_option);
            % else
            srf_hdl=surface(axe_hdl,pnts(:,:,1)',pnts(:,:,2)',pnts(:,:,3)',srf_option);
            zlabel('\itZ');
            % end
            xlabel('\itX');
            ylabel('\itY');

            if nargout > 0, varargout={srf_hdl};end
        end

        function varargout=displayPole(self,axe_hdl,pole_option)
            % draw surface on axes handle
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

            % load properties to calculate u_pole_vctr
            poles=self.getPoles();
            u_deg=self.u_order;v_deg=self.v_order;u_num=self.u_coef_num;v_num=self.v_coef_num;
            u_knotvrtc=self.u_knotvctr;v_knotvrtc=self.v_knotvctr;
            if isempty(u_knotvrtc)
                us=linspace(0,1,u_num)';
            else
                us=interp1(linspace(0,1,u_num-u_deg+1),u_knotvrtc(u_deg+1:u_num+1),linspace(0,1,u_num))';
            end
            u_class=us;if self.sym_x, u_class=us/2+0.5;end
            if isempty(v_knotvrtc)
                vs=linspace(0,1,v_num)';
            else
                vs=interp1(linspace(0,1,v_num-v_deg+1),v_knotvrtc(v_deg+1:v_num+1),linspace(0,1,v_num))';
            end
            v_class=vs;if self.sym_y, v_class=vs/2+0.5;end

            [us,vs]=ndgrid(us,vs);
            [u_class,v_class]=ndgrid(u_class,v_class);

            % calculate class
            cls_x=baseFcnClass(v_class,self.N1x,self.N2x);
            cls_y=baseFcnClass(u_class,self.N1y,self.N2y);
            cls_zv=baseFcnClass(v_class,self.N1zv,self.N2zv);
            cls_zu=baseFcnClass(u_class,self.N1zu,self.N2zu);
            cls=cat(3,cls_x,cls_y,cls_zu.*cls_zv,ones([size(cls_x),self.coef_dim-4]));

            % calculate bias
            if isempty(self.bias_fcn)
                bias=self.B11*(1-us).*(1-vs)+self.B21*us.*(1-vs)+self.B12*(1-us).*vs+self.B22*us.*vs;
                bias=cat(3,zeros([size(bias),2]),bias,zeros([size(bias),self.coef_dim-4]));
            else
                bias=self.bias_fcn(us);
            end

            poles=cls.*poles+bias;
            poles=self.convertLocalToGlobal(poles);

            % draw poles on axe_hdl
            % if self.coef_dim-1 == 2
            %     srf_hdl=surface(axe_hdl,pnts(:,:,1)',pnts(:,:,2)',pole_option);
            % else
            srf_hdl=surface(axe_hdl,poles(:,:,1)',poles(:,:,2)',poles(:,:,3)',pole_option);
            zlabel('\itZ');
            % end
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

            origin=self.origin;
            coord=self.coord;

            % scale of quiver
            scale=norm(diff([axe_hdl.XLim;axe_hdl.YLim;axe_hdl.ZLim],1,2));
            coord=coord*scale*0.05;

            % draw coord on axe_hdl
            % if self.coef_dim-1 == 2
            %     hold on;
            %     quv_hdl(1)=quiver(axe_hdl,origin(1),origin(2),coord(1,1),coord(2,1),0,'Color','r');
            %     quv_hdl(2)=quiver(axe_hdl,origin(1),origin(2),coord(1,2),coord(2,2),0,'Color','g');
            %     hold off;
            % elseif self.coef_dim-1 == 3
            hold on;
            quv_hdl(1)=quiver3(axe_hdl,origin(1),origin(2),origin(3),coord(1,1),coord(2,1),coord(3,1),0,'Color','r');
            quv_hdl(2)=quiver3(axe_hdl,origin(1),origin(2),origin(3),coord(1,2),coord(2,2),coord(3,2),0,'Color','g');
            hold off;
            zlabel('\itZ');
            % end
            xlabel('\itX');
            ylabel('\itY');

            if nargout > 0, varargout={quv_hdl};end
        end
    end

    methods % control surface
        function self=translate(self,tran_vctr)
            % translate curve
            %
            tran_vctr=reshape(tran_vctr,1,[]);
            self.origin=self.origin+tran_vctr(1:self.coef_dim-1);
        end

        function self=rotate(self,rot_mat,rot_cntr)
            % rotate curve
            %
            if nargin < 3 || isempty(rot_cntr), rot_cntr=zeros(1,self.coef_dim-1);end
            rot_cntr=reshape(rot_cntr,1,[]);
            del_vctr=(self.origin-rot_cntr(1:self.coef_dim-1));
            self.origin=rot_cntr(1:self.coef_dim-1)+del_vctr*rot_mat';

            self.coord=rot_mat*self.coord;
        end

        function pnts=convertLocalToGlobal(self,pnts)
            % convert point from local to global coordination
            %
            pnt_size=size(pnts);
            pnts=reshape(pnts,[],self.coef_dim-1);

            % rotate curve
            rot_mat=self.coord;
            pnts=pnts*rot_mat';

            % translate curve
            pnts=pnts+self.origin;

            pnts=reshape(pnts,pnt_size);
        end

        function pnts=convertGlobalToLocal(self,pnts)
            % convert point from global to local coordination
            %
            pnt_size=size(pnts);
            pnts=reshape(pnts,[],self.coef_dim-1);

            % re-translate curve
            pnts=pnts-self.origin;

            % re-rotate curve
            rot_mat=self.coord';
            pnts=pnts*rot_mat';

            pnts=reshape(pnts,pnt_size);
        end

        function srf=convertSpline(self)
            % convert CST Surface into Spline Surface
            %

            % warning('SurfaceCST.convertSpline: can not convert to spline surface exactly');
            u_param=linspace(0,1,21);v_param=linspace(0,1,21);
            msh_1zu=max(min(1/self.N1zu,2),0.5);msh_2zu=max(min(1/self.N2zu,2),0.5);
            if (self.N1zu ~= 0) && (self.N2zu ~= 0) && (~self.sym_x),u_param=u_param.^msh_1zu.*(1-u_param)+u_param.*(1-(1-u_param).^msh_2zu);
            elseif (self.N1zu == 0) && (self.N2zu ~= 0) || (self.sym_x),u_param=1-(1-u_param).^(msh_2zu);
            elseif (self.N1zu ~= 0) && (self.N2zu == 0),u_param=u_param.^(msh_1zu);
            end
            msh_1zv=max(min(1/self.N1zv,2),0.5);msh_2zv=max(min(1/self.N2zv,2),0.5);
            if (self.N1zv ~= 0) && (self.N2zv ~= 0) && (~self.sym_y),v_param=v_param.^msh_1zv.*(1-v_param)+v_param.*(1-(1-v_param).^msh_2zv);
            elseif (self.N1zv == 0) && (self.N2zv ~= 0) || (self.sym_y),v_param=1-(1-v_param).^(msh_2zv);
            elseif (self.N1zv ~= 0) && (self.N2zv == 0),v_param=v_param.^(msh_1zv);
            end
            nodes=self.calPoint({u_param,v_param});
            u_degree=3;v_degree=3;
            srf=interpPointToSurface(nodes,u_degree,v_degree);
        end
    end

    methods % fit surface
        function self=fitSpline(self,nodes,u_degree,v_degree,u_pole_num,v_pole_num,u_nodes,v_nodes)
            % fit Spline base on CST class function
            %
            if nargin < 8
                v_nodes=[];
                if nargin < 7
                    u_nodes=[];
                    if nargin < 6
                        v_pole_num=[];
                        if nargin < 5
                            u_pole_num=[];
                            if nargin < 4
                                v_degree=[];
                                if nargin < 3
                                    u_degree=[];
                                end
                            end
                        end
                    end
                end
            end

            [v_node_num,u_node_num,coef_dim]=size(nodes);
            if coef_dim < 2
                error('SurfaceCST.fitSpline: dimension of nodes can not less than 3')
            end
            if isempty(u_pole_num),u_pole_num=u_node_num;end
            if isempty(v_pole_num),v_pole_num=v_node_num;end
            if u_pole_num > u_node_num || v_pole_num > v_node_num
                error('SurfaceCST.fitSpline: pole_num more than node_num')
            end

            % default value of u_nodes
            if isempty(u_nodes)
                u_nodes=vecnorm(nodes(2:end,:,:)-nodes(1:end-1,:,:),2,3);
                u_nodes=mean(u_nodes,2);u_nodes=[0;cumsum(u_nodes)];
            end
            if all(size(u_nodes) ~= 1), error('SurfaceCST.fitSpline: u_nodes can not be matrix');end
            u_nodes=(u_nodes(:)-min(u_nodes))/(max(u_nodes)-min(u_nodes));
            u_node_knots=interp1(linspace(0,1,length(u_nodes)),u_nodes,linspace(0,1,u_pole_num));

            % default value of v_nodes
            if isempty(v_nodes)
                v_nodes=vecnorm(nodes(:,2:end,:)-nodes(:,1:end-1,:),2,3);
                v_nodes=mean(v_nodes,1);v_nodes=[0;cumsum(v_nodes')];
            end
            if all(size(v_nodes) ~= 1), error('SurfaceCST.fitSpline: v_nodes can not be matrix');end
            v_nodes=(v_nodes(:)-min(v_nodes))/(max(v_nodes)-min(v_nodes));
            v_node_knots=interp1(linspace(0,1,length(v_nodes)),v_nodes,linspace(0,1,v_pole_num));

            u_mults=[u_degree+1,ones(1,u_pole_num-u_degree-1),u_degree+1];
            u_knots=linspace(0,1,u_pole_num-u_degree+1);
            for j=2:u_pole_num-u_degree
                u_knots(j)=mean(u_node_knots(j:j+u_degree-1));
            end
            % modify
            u_knots=interp1(linspace(0,1,length(u_knots)),u_knots,linspace(0,1,u_pole_num-u_degree+1));
            u_list=baseKnotVctr(u_mults,u_knots);

            v_mults=[v_degree+1,ones(1,v_pole_num-v_degree-1),v_degree+1];
            v_knots=linspace(0,1,v_pole_num-v_degree+1);
            for j=2:v_pole_num-v_degree
                v_knots(j)=mean(v_node_knots(j:j+v_degree-1));
            end
            % modify
            v_knots=interp1(linspace(0,1,length(v_knots)),v_knots,linspace(0,1,v_pole_num-v_degree+1));
            v_list=baseKnotVctr(v_mults,v_knots);

            % translate node to local coordinate
            nodes=self.convertGlobalToLocal(nodes);

            % remove bias
            dim_class=3;
            if isempty(self.bias_fcn)
                bias=self.B11*(1-u_nodes).*(1-v_nodes)+self.B21*u_nodes.*(1-v_nodes)+self.B12*(1-u_nodes).*v_nodes+self.B22*u_nodes.*v_nodes;
            else
                bias=self.bias_fcn([u_nodes,v_nodes]);
            end
            nodes(:,:,dim_class)=nodes(:,:,dim_class)-bias;

            % base on node point list inverse calculate control point list
            u_fit_mat=zeros(u_node_num,u_pole_num);
            [N_list,idx_srt,idx_end]=baseFcnN(u_nodes,u_degree,u_list);
            for deg_idx=1:u_degree+1
                idx=sub2ind([u_node_num,u_pole_num],(1:(u_node_num))',idx_srt+(deg_idx-1));
                u_fit_mat(idx)=N_list(:,deg_idx);
            end
            v_fit_mat=zeros(v_node_num,v_pole_num);
            [N_list,idx_srt,idx_end]=baseFcnN(v_nodes,v_degree,v_list);
            for deg_idx=1:v_degree+1
                idx=sub2ind([v_node_num,v_pole_num],(1:(v_node_num))',idx_srt+(deg_idx-1));
                v_fit_mat(idx)=N_list(:,deg_idx);
            end

            % add class coefficient
            v_class=v_nodes;
            if self.sym_y,v_class=(v_nodes/2)+0.5;end
            u_class=u_nodes;
            if self.sym_x,u_class=(u_nodes/2)+0.5;end

            v_fit_mat_class_x=v_fit_mat.*baseFcnClass(v_class,self.N1x,self.N2x);
            v_fit_mat_class_y=v_fit_mat;
            v_fit_mat_class_zv=v_fit_mat.*baseFcnClass(v_class,self.N1zv,self.N2zv);

            u_fit_mat_class_x=u_fit_mat;
            u_fit_mat_class_y=u_fit_mat.*baseFcnClass(u_class,self.N1y,self.N2y);
            u_fit_mat_class_zu=u_fit_mat.*baseFcnClass(u_class,self.N1zu,self.N2zu);

            v_fit_mat=v_fit_mat';
            v_fit_mat_class_x=v_fit_mat_class_x';
            v_fit_mat_class_y=v_fit_mat_class_y';
            v_fit_mat_class_zv=v_fit_mat_class_zv';

            % reverse calculate control point
            poles(:,:,1)=u_fit_mat_class_x\nodes(:,:,1)/v_fit_mat_class_x;
            poles(:,:,2)=u_fit_mat_class_y\nodes(:,:,2)/v_fit_mat_class_y;
            poles(:,:,3)=u_fit_mat_class_zu\nodes(:,:,3)/v_fit_mat_class_zv;
            if coef_dim > 3
                poles(:,:,4:end)=u_fit_mat\nodes(:,:,4:end)/v_fit_mat;
            end
            weights=ones(v_pole_num,u_pole_num);

            % add spline and shape function
            self=self.addSpline(poles,u_degree,v_degree,u_mults,v_mults,u_knots,v_knots,weights);

            % fit data
            self.nodes=nodes;
            self.u_nodes=u_nodes;
            self.v_nodes=v_nodes;
            self.u_fit_mat=u_fit_mat;
            self.v_fit_mat=v_fit_mat;
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
            uv_min=[0+geom_tol,0+geom_tol];uv_max=[1-geom_tol,1-geom_tol];

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
            u_base_list=base_list;
            v_base_list=base_list;
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

%% function

function [pnts,wts]=calSurface(coef_dim,u_coef_num,v_coef_num,coefs,u_knots,v_knots,u_order,v_order,us,vs)
% calculate spline point and weight
%
if nargin < 10,vs=[];end
if isempty(vs),uvs=us;end

if iscell(us)
    us=uvs{1}(:);vs=uvs{2}(:);
else
    uvs=cat(3,us,vs);
    siz=size(uvs);uvs=reshape(uvs,[],2);
    us=uvs(:,1);vs=uvs(:,2);
end

if iscell(uvs)
    % input is {u_list,v_list}, mesh point
    u_num=length(us);v_num=length(vs);

    % evaluate along the u direction
    coefs_u=reshape(coefs,[u_coef_num,v_coef_num*coef_dim]);
    if u_coef_num == 2
        pnts=coefs_u(1,:).*(1-us)+coefs_u(2,:).*us;
    else
        min_u=u_knots(1);max_u=u_knots(end);
        us_k=us*(max_u-min_u)+min_u;
        [u_N_list,u_idx_srt,~]=baseFcnN(us_k,u_order,u_knots);
        pnts=zeros(length(us_k),size(coefs_u,2));
        for deg_idx=1:u_order+1
            pnts=pnts+u_N_list(:,deg_idx).*coefs_u(u_idx_srt+(deg_idx-1),:);
        end
    end
    pnts=reshape(pnts,[u_num,v_coef_num,coef_dim]);

    % evaluate along the v direction
    coefs_v=permute(pnts,[2,1,3]);
    coefs_v=reshape(coefs_v,v_coef_num,coef_dim*u_num);
    if v_coef_num == 2
        pnts=coefs_v(1,:).*(1-vs)+coefs_v(2,:).*vs;
    else
        min_v=v_knots(1);max_v=v_knots(end);
        vs_k=vs*(max_v-min_v)+min_v;
        [v_N_list,v_idx_srt,~]=baseFcnN(vs_k,v_order,v_knots);
        pnts=zeros(length(vs_k),size(coefs_v,2));
        for deg_idx=1:v_order+1
            pnts=pnts+v_N_list(:,deg_idx).*coefs_v(v_idx_srt+(deg_idx-1),:);
        end
    end
    pnts=reshape(pnts,[v_num,u_num,coef_dim]);
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
    coefs_u=reshape(coefs,[u_coef_num,v_coef_num*coef_dim]);
    if u_coef_num == 2
        pnts=coefs_u(1,:).*(1-us)+coefs_u(2,:).*us;
    else
        min_u=u_knots(1);max_u=u_knots(end);
        us_k=us*(max_u-min_u)+min_u;
        [u_N_list,u_idx_srt,~]=baseFcnN(us_k,u_order,u_knots);
        pnts=zeros(length(us_k),size(coefs_u,2));
        for deg_idx=1:u_order+1
            pnts=pnts+u_N_list(:,deg_idx).*coefs_u(u_idx_srt+(deg_idx-1),:);
        end
    end
    pnts=reshape(pnts,[pnt_num,v_coef_num,coef_dim]);

    % evaluate along the v direction
    coefs_v=reshape(pnts,pnt_num*v_coef_num,coef_dim);
    if v_coef_num == 2
        pnts=coefs_v(1:pnt_num,:).*(1-vs)+coefs_v(pnt_num+1:end,:).*vs;
    else
        min_v=v_knots(1);max_v=v_knots(end);
        vs_k=vs*(max_v-min_v)+min_v;
        [v_N_list,v_idx_srt,~]=baseFcnN(vs_k,v_order,v_knots);
        pnts=zeros(pnt_num,coef_dim);
        for deg_idx=1:v_order+1
            idx=sub2ind([pnt_num,v_coef_num],(1:pnt_num)',v_idx_srt+(deg_idx-1));
            pnts=pnts+v_N_list(:,deg_idx).*coefs_v(idx,:);
        end
    end

    wts=pnts(:,end);
    pnts=pnts(:,1:end-1);
    if nargout < 2
        pnts=pnts./wts;
    end

    pnts=reshape(pnts,[siz(1:end-1),coef_dim-1]);
    wts=reshape(wts,[siz(1:end-1),1]);
end
end

function [pnts,dpnts_du,dpnts_dv]=differFcn2Robust(fcn,us,vs,step,dim)
% calculate function gradient by robust differ
%
pnts=reshape(fcn(us,vs),[],dim);

pnts_uf=reshape(fcn(min(us+step,1),vs),[],dim);
pnts_ub=reshape(fcn(max(us-step,0),vs),[],dim);
bools_f=(us+step) >= 1;
bools_b=(us-step) <= 0;
bools_c=~any([bools_f,bools_b],2);
dpnts_du=zeros(size(pnts)); % allocate memory
dpnts_du(bools_c,:)=(pnts_uf(bools_c,:)-pnts_ub(bools_c,:))/2/step;
dpnts_du(bools_f,:)=(pnts(bools_f,:)-pnts_ub(bools_f,:))/step;
dpnts_du(bools_b,:)=(pnts_uf(bools_b,:)-pnts(bools_b,:))/step;
dpnts_du=real(dpnts_du);

pnts_vf=reshape(fcn(us,min(vs+step,1)),[],dim);
pnts_vb=reshape(fcn(us,max(vs-step,0)),[],dim);
bools_f=(vs+step) >= 1;
bools_b=(vs-step) <= 0;
bools_c=~any([bools_f,bools_b],2);
dpnts_dv=zeros(size(pnts)); % allocate memory
dpnts_dv(bools_c,:)=(pnts_vf(bools_c,:)-pnts_vb(bools_c,:))/2/step;
dpnts_dv(bools_f,:)=(pnts(bools_f,:)-pnts_vb(bools_f,:))/step;
dpnts_dv(bools_b,:)=(pnts_vf(bools_b,:)-pnts(bools_b,:))/step;
dpnts_dv=real(dpnts_dv);
end
