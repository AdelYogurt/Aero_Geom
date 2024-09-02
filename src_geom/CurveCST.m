classdef CurveCST
    % Class-Shape Transformation (CST) Curve
    %
    % notice:
    % origin reference of shape function is local coordinate
    % point calculate by class and shape will be translate to global
    %
    % reference:
    % [1] Kulfan B. A Universal Parametric Geometry Representation Method -
    % "CST" [C]. 45th AIAA Aerospace Sciences Meeting and Exhibit, Reno,
    % Nevada, U.S.A.
    %
    properties % CST parameter
        form='.CST.'; % curve_form

        % class parameter
        N1=[];
        N2=[];

        % shape parameter
        coef_dim=3; % curve dimension of coefs
        u_coef_num=2; % curve number of coefs
        coefs=[0,1,1;1,1,1]; % curve coefs matrix, equal to permute(cat(2,poles.*weights,weights),[2,1]);
        u_knotvctr=[0,0,1,1]; % curve knot_vector
        u_order=1; % curve degree

        shape_fcn=[]; % function handle

        % bias parameter
        B1=0;
        B2=0;

        bias_fcn=[]; % function handle

        coord=[]; % coordination of curve, [Ax1,Ax2]
        origin=[]; % origin point of coordination
        sym=false; % if true, u_class will start from 0.5 to 1

        deriv_crv=[];
    end

    properties(Access=private) % fit point set parameter
        nodes; % fit point set
        u_nodes; % fit point parameter
        fit_mat; % fit matrix
    end

    methods % define curve
        function self=CurveCST(C_par,S_par,B_par,sym)
            % generate CST curve
            %
            % input:
            % C_par (list): N1 and N2 array([N1, N2]), default is [0.0,0.0]
            % S_par (list): LX and LY array([LX, LY]), default is [1.0,1.0]
            % B_par (list): B1 and B2 array([B1, B2]), default is [0.0,0.0]
            % sym (boolean): whether using half of u while calculate class_fcn
            %
            % output:
            % CurveCST
            %
            % notice:
            % if C_par is [0.0, 0.0] or empty, class_fcn will equal to 1.0
            %
            if nargin < 4
                sym=[];
                if nargin < 3
                    B_par=[];
                    if nargin < 2
                        S_par=[];
                        if nargin < 1
                            C_par=[];
                        end
                    end
                end
            end

            % default value
            if isempty(sym),sym=false;end
            if isempty(B_par),B_par=[0.0,0.0];end
            if isempty(S_par),S_par=[1.0,1.0];end
            if isempty(C_par),C_par=[0.0,0.0];end

            self.sym=sym;

            % process C_par
            self.N1=C_par(1);
            self.N2=C_par(2);

            % process S_par
            if isa(S_par,'function_handle')
                shape_fcn=S_par;
                self.shape_fcn=shape_fcn;
                pnt=shape_fcn(0);
                coef_dim=length(pnt)+1;

                self.coef_dim=coef_dim;
            else
                LX=S_par(1);LY=S_par(2);

                pnt_1=[0,LY,1];
                pnt_2=[LX,LY,1];
                coef_dim=3;number=2;order=1;
                coefs=[pnt_1;pnt_2];
                u_knots=[0,0,1,1];

                self.coef_dim=coef_dim;
                self.u_coef_num=number;
                self.coefs=coefs;
                self.u_order=order;
                self.u_knotvctr=u_knots;
            end

            % process B_par
            if isa(B_par,'function_handle')
                self.bias_fcn=B_par;
            else
                self.B1=B_par(1);
                self.B2=B_par(2);
            end

            self.coord=eye(coef_dim-1);
            self.origin=zeros(1,coef_dim-1);
        end
        
        function self=deriv(self,~)
            % generate derivation data
            %
            if isempty(self.shape_fcn)
                if isempty(self.deriv_crv)
                    if self.u_coef_num == 2
                        dcoefs_du=self.coefs(2,:)-self.coefs(1,:);
                        dorder_du=0;dknots_du=[0,1];
                    elseif self.u_order == (self.u_coef_num-1) && isempty(self.u_knotvctr)
                        dcoefs_du=self.u_order*diff(self.coefs,1,1);
                        dorder_du=self.u_order-1;dknots_du=self.u_knotvctr(2:end-1);
                    elseif ~isempty(self.u_knotvctr)
                        [dcoefs_du,dorder_du,dknots_du]=BSpline.deriv(self.coefs,self.u_order,self.u_knotvctr);
                    else
                        error('CurveCST.deriv: error shape parameter define')
                    end

                    self.deriv_crv.coef_dim=self.coef_dim;
                    self.deriv_crv.u_coef_num=self.u_coef_num-1;
                    self.deriv_crv.coefs=dcoefs_du;
                    self.deriv_crv.u_knots=dknots_du;
                    self.deriv_crv.u_order=dorder_du;
                end
            end
        end

        function self=addSpline(self,poles,u_degree,u_mults,u_knots,weights)
            % add Bezier/B-Splines Curve as shape function
            %
            % calling:
            % crv=crv.addSpline(poles,degree,mults,knots,weights)
            % crv=crv.addSpline(coefs,knots)
            %
            % input:
            % poles (matrix): control point, pole_num x dimension matrix
            % degree (matrix): optional input
            % mults (matrix): optional input
            % knots (matrix): optional input
            % weights (matrix): optional input
            %
            % output:
            % CurveCST
            %
            % notice:
            % degree default is pole_num-1 which will be Bezier curve
            %
            if nargin < 6
                weights=[];
                if nargin < 5
                    u_knots=[];
                    if nargin < 4
                        u_mults=[];
                        if nargin < 3
                            u_degree=[];
                        end
                    end
                end
            end

            [u_pole_num,coef_dim]=size(poles);

            % default value
            if isempty(u_degree), u_degree=u_pole_num-1;end

            % different input mode
            if length(u_degree) > 1
                coefs=poles;
                u_knotvrtc=u_degree;
                u_knotvrtc=sort(u_knotvrtc);
                u_order=length(u_knotvrtc)-u_pole_num-1;
            else
                if u_pole_num < (u_degree+1)
                    error('CurveCST.addSpline: pole_num less than degree+1');
                end

                % u default value
                if isempty(u_mults) && isempty(u_knots)
                    u_mults=[u_degree+1,ones(1,u_pole_num-u_degree-1),u_degree+1];
                    u_knots=linspace(0,1,u_pole_num-u_degree+1);
                elseif ~isempty(u_mults) && isempty(u_knots)
                    u_knots=linspace(0,1,length(u_mults));
                elseif isempty(u_mults) && ~isempty(u_knots)
                    error('CurveCST.addSpline: need mults input');
                end
                u_knotvrtc=baseKnotVctr(u_mults,u_knots);
                u_knotvrtc=sort(u_knotvrtc);

                if length(u_knotvrtc) ~= u_pole_num+u_degree+1
                    error('CurveCST.addSpline: knot_num is not equal to pole_num+degree+1');
                end

                if min(u_knotvrtc) ~= 0.0 || max(u_knotvrtc) ~= 1.0
                    warning('CurveCST.addSpline: knots will be rescale to [0.0, 1.0]')
                    u_knotvrtc=(u_knotvrtc-min(u_knotvrtc))/(max(u_knotvrtc)-min(u_knotvrtc));
                end

                if isempty(weights), weights=ones(1,u_pole_num);end
                if ~isempty(weights), weights=weights(:);end

                u_order=u_degree;
                coefs=cat(2,poles.*weights,weights);
                coef_dim=coef_dim+1;
            end

            self.coef_dim=coef_dim;
            self.u_coef_num=u_pole_num;
            self.u_order=u_order;
            self.coefs=coefs;
            self.u_knotvctr=u_knotvrtc;
        end
        
        function [poles,weights]=getPoles(self)
            % convert coefs to poles and weights
            %
            poles=self.coefs(:,1:end-1);
            weights=self.coefs(:,end);
            weights(weights == 0)=1;
            poles=poles./weights;
        end

        function [pnts,cls,sps,bias]=calPoint(self,us)
            % calculate point on curve
            %

            % calculate shape
            if isempty(self.shape_fcn)
                sps=calCurve(self.coef_dim,self.u_coef_num,self.coefs,self.u_knotvctr,self.u_order,us);
            else
                sps=self.shape_fcn(us);
            end
            
            % calculate perpare
            us=us(:);
            u_class=us;if self.sym,u_class=(u_class/2)+0.5;end

            % calculate class
            cls=baseFcnClass(u_class,self.N1,self.N2);
            cls=[ones(size(cls)),cls,ones(size(cls,1),self.coef_dim-3)];

            % calculate bias
            if isempty(self.bias_fcn)
                bias=self.B1*(1-us)+self.B2*us;
                bias=[zeros(size(bias)),bias,zeros(size(bias,1),self.coef_dim-3)];
            else
                bias=self.bias_fcn(us);
                bias=[zeros(size(bias)),bias,zeros(size(bias,1),self.coef_dim-3)];
            end

            % calculate point
            pnts=cls.*sps+bias;
            pnts=self.convertLocalToGlobal(pnts);
        end

        function [pnts,dpnts_du]=calGradient(self,us,step)
            % calculate gradient of curve
            %
            if nargin < 3 || isempty(step),step=sqrt(eps);end
            self=self.deriv(1);

            % calculate shape
            if isempty(self.shape_fcn)
                [sps_p,sps_w]=calCurve(self.coef_dim,self.u_coef_num,self.coefs,self.u_knotvctr,self.u_order,us);
                dcrv_du=self.deriv_crv;
                [dsps_du_p,dsps_du_w]=calCurve(dcrv_du.coef_dim,dcrv_du.u_coef_num,dcrv_du.coefs,dcrv_du.u_knots,dcrv_du.u_order,us);
                sps=sps_p./sps_w;
                dsps_du=(dsps_du_p-sps.*dsps_du_w)./sps_w;

                du=self.u_knotvctr(end)-self.u_knotvctr(1);
                if self.u_coef_num ~= 2,dsps_du=dsps_du*du;end
            else
                [sps,dsps_du]=differFcn1Robust(self.shape_fcn,us,step);
            end

            % calculate perpare
            us=us(:);
            u_class=us;if self.sym,u_class=(u_class/2)+0.5;end

            % calculate class
            [cls,dcls_du]=baseFcnClass(u_class,self.N1,self.N2);
            if self.sym,dcls_du=dcls_du/2;end,dcls_du=min(dcls_du,1e6);
            cls=[ones(size(cls)),cls,ones(size(cls,1),self.coef_dim-3)];
            dcls_du=[zeros(size(dcls_du)),dcls_du,zeros(size(dcls_du,1),self.coef_dim-3)];

            % calculate bias
            if isempty(self.bias_fcn)
                bias=self.B1*(1-us)+self.B2*us;
                dbias_du=self.B2-self.B1;
                bias=[zeros(size(bias)),bias,zeros(size(bias,1),self.coef_dim-3)];
                dbias_du=[zeros(size(dbias_du)),dbias_du,zeros(size(dbias_du,1),self.coef_dim-3)];
            else
                [bias,dbias_du]=differFcn1Robust(self.bias_fcn,us,step);
            end

            % calculate gradient
            pnts=cls.*sps+bias;
            dpnts_du=dcls_du.*sps+cls.*dsps_du+dbias_du;
            pnts=self.convertLocalToGlobal(pnts);
            dpnts_du=dpnts_du*self.coord';
        end

        function [k1,k2,x1,x2]=calTangTorl(self,torl)
            % calculate differ tangent by control discrete error in torl
            %
            CN=(self.N1./(self.N1+self.N2)).^self.N1.*(self.N2./(self.N1+self.N2)).^self.N2;
            if isempty(self.shape_fcn)
                poles=self.getPoles();
                K1=poles(1,2);
                K2=poles(end,2);
            else
                K1=self.shape_fcn(0);K1=K1(2);
                K2=self.shape_fcn(1);K2=K2(2);
            end
            KY1=K1/CN;
            KY2=K2/CN;
            KX=abs(K2(1)-K1(1));
            if self.sym,KX=KX*2;end

            if self.N1 == 1
                k1=self.N1*KY1/KX;
                x1=0;
            else
                k1=self.N1*KY1/KX*(torl/abs(KY1)/abs(1-self.N1))^((self.N1-1)/self.N1);
                x1=KX*(k1*KX/self.N1/abs(KY1))^(1/(self.N1-1));
            end

            if self.N2 == 1
                k2=self.N2*KY2/KX;
                x2=0;
            else
                k2=self.N2*KY2/KX*(torl/abs(KY2)/abs(1-self.N2))^((self.N2-1)/self.N2);
                x2=KX*(k1*KX/self.N2/abs(KY2))^(1/(self.N2-1));
            end
        end

        function varargout=displayGeom(self,axe_hdl,crv_option,param)
            % display curve on axes
            %
            if nargin < 4
                param=[];
                if nargin < 3
                    crv_option=[];
                    if nargin < 2
                        axe_hdl=[];
                    end
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end
            if isempty(crv_option), crv_option=struct();end
            if isempty(param),param=51;end
            if length(param) == 1,param=linspace(0,1,param);end
            msh_1=max(min(1/self.N1,2),0.5);msh_2=max(min(1/self.N2,2),0.5);
            if (self.N1 ~= 0) && (self.N2 ~= 0) && (~self.sym),param=param.^msh_1.*(1-param)+param.*(1-(1-param).^msh_2);
            elseif (self.N1 == 0) && (self.N2 ~= 0) || (self.sym),param=1-(1-param).^(msh_2);
            elseif (self.N1 ~= 0) && (self.N2 == 0),param=param.^(msh_1);
            end
            pnts=self.calPoint(param);

            % draw pnts on axe_hdl
            if self.coef_dim-1 == 2
                ln_hdl=line(axe_hdl,pnts(:,1),pnts(:,2),crv_option);
            elseif self.coef_dim-1 == 3
                ln_hdl=line(axe_hdl,pnts(:,1),pnts(:,2),pnts(:,3),crv_option);
                zlabel('z');
            end
            xlabel('x');
            ylabel('y');

            if nargout > 0, varargout={ln_hdl};end
        end

        function varargout=displayPoles(self,axe_hdl,pole_option)
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
                pole_option=struct('Marker','s','LineStyle','--','Color','r');
            end

            % load properties to calculate u_pole_vctr
            poles=self.getPoles();
            deg=self.u_order;num=self.u_coef_num;
            knotvrtc=self.u_knotvctr;
            min_u=self.u_knotvctr(1);max_u=self.u_knotvctr(end);du=max_u-min_u;
            us=interp1(linspace(0,1,num-deg+1),knotvrtc(deg+1:num+1),linspace(min_u,max_u,num))';
            us=us(:);us_norm=(us-min_u)/du;u_class=us_norm;
            if self.sym,u_class=(u_class/2)+0.5;end

            % calculate class
            cls=baseFcnClass(u_class,self.N1,self.N2);
            cls=[ones(size(cls)),cls,ones(size(cls,1),self.coef_dim-3)];

            % calculate bias
            if isempty(self.bias_fcn)
                bias=self.B1*(1-us_norm)+self.B2*us_norm;
                bias=[zeros(size(bias)),bias,zeros(size(bias,1),self.coef_dim-3)];
            else
                bias=self.bias_fcn(us);
            end

            poles=cls.*poles+bias;
            poles=self.convertLocalToGlobal(poles);

            % draw poles on axe_hdl
            if self.coef_dim-1 == 2
                ln_hdl=line(axe_hdl,poles(:,1),poles(:,2),pole_option);
            elseif self.coef_dim-1 == 3
                ln_hdl=line(axe_hdl,poles(:,1),poles(:,2),poles(:,3),pole_option);
                zlabel('z');
            end
            xlabel('x');
            ylabel('y');

            if nargout > 0, varargout={ln_hdl};end
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
            if self.coef_dim-1 == 2
                hold on;
                quv_hdl(1)=quiver(axe_hdl,origin(1),origin(2),coord(1,1),coord(2,1),0,'Color','r');
                hold off;
            elseif self.coef_dim-1 == 3
                hold on;
                quv_hdl(1)=quiver3(axe_hdl,origin(1),origin(2),origin(3),coord(1,1),coord(2,1),coord(3,1),0,'Color','r');
                hold off;
                zlabel('z');
            end
            xlabel('x');
            ylabel('y');

            if nargout > 0, varargout={quv_hdl};end
        end
    end

    methods % control curve
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

        function self=insertDimension(self,dim_add,bias)
            % add additional geometry dimension
            %
            if nargin < 3 || isempty(bias), bias=0;end
            self.coefs=[self.coefs(:,1:dim_add),...
                ones(self.u_coef_num,1)*bias,...
                self.coefs(:,(dim_add+1):end)];
            self.coef_dim=self.coef_dim+1;
        end

        function pnts=convertLocalToGlobal(self,pnts)
            % convert point from local to global coordination
            %

            % rotate curve
            rot_mat=self.coord;
            pnts=pnts*rot_mat';

            % translate curve
            pnts=pnts+self.origin;
        end

        function pnts=convertGlobalToLocal(self,pnts)
            % convert point from global to local coordination
            %

            % re-translate curve
            pnts=pnts-self.origin;

            % re-rotate curve
            rot_mat=self.coord';
            pnts=pnts*rot_mat';
        end

        function crv=convertSpline(self)
            % convert CST Curve to Spline Curve
            %
            % reference: [1] Marshall D D. Creating Exact Bezier
            % Representations of CST Shapes [C]. 21st AIAA Computational
            % Fluid Dynamics Conference, San Diego, CA.
            %

            min_u=self.u_knotvctr(1);max_u=self.u_knotvctr(end);
            if isempty(self.shape_fcn) && isempty(self.bias_fcn) && self.u_order == self.u_coef_num-1
                % convert CST into Bezier excatly
                CN=(self.N1./(self.N1+self.N2)).^self.N1.*(self.N2./(self.N1+self.N2)).^self.N2;
                CN((self.N1 == 0) & (self.N2 == 0))=1;

                if (self.N1 == 0.5 && self.N2 == 1.0) || ...
                    (self.N1 == 1.0 && self.N2 == 1.0) || ...
                    (self.N1 == 1.0 && self.N2 == 0.0) || ...
                    (self.N1 == 0.0 && self.N2 == 0.0) || ...
                    (self.N1 == 0.0 && self.N2 == 1.0)

                    poles=self.getPoles();
                    deg=size(poles,1)-1;
                    % calculate a list
                    a_list=convertBezierToPoly(poles);
                    if self.N1 == 0.5 && self.N2 == 1.0
                        % calculate b list
                        b_list=zeros(2*deg+3+1,size(poles,2));

                        % other dimension
                        dim_idx=1:size(poles,2);dim_idx(2)=[];
                        i=0:2:2*deg;
                        b_list(i+1,dim_idx)=a_list(i/2+1,dim_idx);

                        % class dimension
                        b_list(0+1,2)=self.B1;
                        b_list(1+1,2)=a_list(0+1,2)/CN;
                        b_list(2+1,2)=self.B2-self.B1;
                        i=3:2:2*deg+1;
                        b_list(i+1,2)=(a_list((i-1)/2+1,2)-a_list((i-1)/2,2))/CN;
                        b_list(2*deg+3+1,2)=-a_list(deg+1,2)/CN;
                    elseif self.N1 == 1.0 && self.N2 == 1.0
                        % calculate b list
                        b_list=zeros(deg+2+1,size(poles,2));

                        % other dimension
                        dim_idx=1:size(poles,2);dim_idx(2)=[];
                        i=0:deg;
                        b_list(i+1,dim_idx)=a_list(i+1,dim_idx);

                        % class dimension
                        b_list(0+1,2)=self.B1;
                        b_list(1+1,2)=a_list(0+1,2)/CN+self.B2-self.B1;
                        i=2:deg+1;
                        b_list(i+1,2)=(a_list(i-1+1,2)-a_list(i-2+1,2))/CN;
                        b_list(deg+2+1,2)=-a_list(deg+1,2)/CN;
                    elseif self.N1 == 1.0 && self.N2 == 0.0
                        % calculate b list
                        b_list=zeros(deg+1+1,size(poles,2));

                        % other dimension
                        dim_idx=1:size(poles,2);dim_idx(2)=[];
                        i=0:deg;
                        b_list(i+1,dim_idx)=a_list(i+1,dim_idx);

                        % class coef_dim
                        b_list(0+1,2)=self.B1;
                        b_list(1+1,2)=a_list(0+1,2)/CN+self.B2-self.B1;
                        i=2:deg+1;
                        b_list(i+1,2)=a_list(i-1+1,2)/CN;
                    elseif self.N1 == 0.0 && self.N2 == 0.0
                        % calculate b list
                        b_list=zeros(deg+1,size(poles,2));

                        % other coef_dim
                        dim_idx=1:size(poles,2);dim_idx(2)=[];
                        i=0:deg;
                        b_list(i+1,dim_idx)=a_list(i+1,dim_idx);

                        % class coef_dim
                        b_list(0+1,2)=a_list(0+1,2)+self.B1;
                        b_list(1+1,2)=a_list(1+1,2)+self.B2-self.B1;
                        i=2:deg;
                        b_list(i+1,2)=a_list(i+1,2);
                    elseif self.N1 == 0.0 && self.N2 == 1.0
                        % calculate b list
                        b_list=zeros(deg+1+1,size(poles,2));

                        % other coef_dim
                        dim_idx=1:size(poles,2);dim_idx(2)=[];
                        i=0:deg;
                        b_list(i+1,dim_idx)=a_list(i+1,dim_idx);

                        % class coef_dim
                        b_list(0+1,2)=a_list(0+1,2)/CN+self.B1;
                        b_list(1+1,2)=(a_list(1+1,2)-a_list(0+1,2))/CN+self.B2-self.B1;
                        i=2:deg;
                        b_list(i+1,2)=(a_list(i+1,2)-a_list(i-1+1,2))/CN;
                        b_list(deg+1+1,2)=-a_list(deg+1,2)/CN;
                    end
                    % calculate pole list
                    poles=convertPolyToBezier(b_list);
                    poles=self.convertLocalToGlobal(poles);
                    crv=Curve(poles);
                    crv.u_knotvctr=crv.u_knotvctr*(max_u-min_u)+min_u;
                    return;
                end
            end

            % warning('CurveCST.convertSpline: can not convert to spline curve exactly');
            param=linspace(0,1,21);
            msh_1=max(min(1/self.N1,2),0.5);msh_2=max(min(1/self.N2,2),0.5);
            if (self.N1 ~= 0) && (self.N2 ~= 0) && (~self.sym),param=param.^msh_1.*(1-param)+param.*(1-(1-param).^msh_2);
            elseif (self.N1 == 0) && (self.N2 ~= 0) || (self.sym),param=1-(1-param).^(msh_2);
            elseif (self.N1 ~= 0) && (self.N2 == 0),param=param.^(msh_1);
            end
            pnts=self.calPoint(param);
            degree=3;
            crv=interpPointToCurve(pnts,degree);
            crv.u_knotvctr=crv.u_knotvctr*(max_u-min_u)+min_u;
        end
    end

    methods % fit curve
        function self=fitSpline(self,nodes,degree,pole_num,u_nodes)
            % fit Spline base on CST class function
            %
            if nargin < 5
                u_nodes=[];
                if nargin < 4
                    pole_num = [];
                    if nargin < 3
                        degree=[];
                    end
                end
            end

            [node_num,coef_dim]=size(nodes);
            if coef_dim < 2
                error('CurveCST.fitSpline: dimension of nodes can not less than 2')
            end
            if isempty(pole_num),pole_num=node_num;end
            if pole_num > node_num
                error('CurveCST.fitSpline: pole_num can not more than node_num')
            end

            % default value of u_nodes vector
            if isempty(u_nodes)
                u_nodes=vecnorm(nodes(2:end,:)-nodes(1:end-1,:),2,2);
                u_nodes=[0;cumsum(u_nodes)];
            end
            u_nodes=(u_nodes(:)-min(u_nodes))/(max(u_nodes)-min(u_nodes));
            u_node_knots=interp1(linspace(0,1,length(u_nodes)),u_nodes,linspace(0,1,pole_num));

            % reference:
            % [1] 施法中. 计算机辅助几何设计与非均匀有理B样条[M]. 208-281.
            % [2] https://blog.csdn.net/he_nan/article/details/107134271
            mults=[degree+1,ones(1,pole_num-degree-1),degree+1];
            knots=linspace(0,1,pole_num-degree+1);
            for j=2:pole_num-degree
                knots(j)=mean(u_node_knots(j:j+degree-1));
            end
            knots=interp1(linspace(0,1,length(knots)),knots,linspace(0,1,pole_num-degree+1)); % modify
            u_list=baseKnotVctr(mults,knots);

            % translate node to local coordinate
            nodes=self.convertGlobalToLocal(nodes);

            % remove bias
            dim_class=2;
            if isempty(self.bias_fcn)
                bias=self.B1*(1-u_nodes)+self.B2*u_nodes;
            else
                bias=self.bias_fcn(u_nodes);
            end
            nodes(:,dim_class)=nodes(:,dim_class)-bias;

            % base on node point list inverse calculate control point list
            fit_mat=zeros(node_num,pole_num);
            [N_list,idx_srt,idx_end]=baseFcnN(u_nodes,degree,u_list);
            for deg_idx=1:degree+1
                idx=sub2ind([node_num,pole_num],(1:(node_num))',idx_srt+(deg_idx-1));
                fit_mat(idx)=N_list(:,deg_idx);
            end

            % add class coefficient
            u_class=u_nodes;
            if self.sym,u_class=(u_nodes/2)+0.5;end
            fit_mat_class=fit_mat.*baseFcnClass(u_class,self.N1,self.N2);

            % reverse calculate control point
            poles=zeros(pole_num,size(nodes,2));
            dim_idx=1:size(nodes,2);dim_idx(dim_class)=[];
            poles(:,dim_idx)=fit_mat\nodes(:,dim_idx);
            poles(:,dim_class)=fit_mat_class\nodes(:,dim_class);

            % add spline and shape function
            self=self.addSpline(poles,degree,mults,knots);

            % fit data
            self.nodes=nodes;
            self.u_nodes=u_nodes;
            self.fit_mat=fit_mat;
        end

        function [fit_err,self]=optimClass(self,optim_option)
            % optimization coefficient of class function
            %
            if nargin < 2
                optim_option=[];
            end

            if isempty(optim_option)
                optim_option=optimoptions('fminunc','Display','none','FiniteDifferenceStepSize',1e-5);
            end

            C_par_low_bou=[0,0];
            C_par_up_bou=[1e3,1e3];
            C_par_init=[self.N1,self.N2];
            obj_fit=@(C_par) self.fitError(C_par);
            [C_par,~]=fminunc(obj_fit,C_par_init,optim_option);

            C_par=max(C_par,C_par_low_bou);
            C_par=min(C_par,C_par_up_bou);
            [fit_err,self]=fitError(self,C_par);
        end

        function [RMSE,self]=fitError(self,C_par)
            % fit error of C_par
            %

            % load fit data
            u_nodes=self.u_nodes;
            nodes=self.nodes;
            fit_matrix=self.fit_mat;

            if nargin > 1
                C_par_low_bou=[0,0];
                C_par_up_bou=[1e3,1e3];
                C_par=max(C_par,C_par_low_bou);
                C_par=min(C_par,C_par_up_bou);

                % updata N1, N2
                self.N1=C_par(1);self.N2=C_par(2);

                pole_num=size(fit_matrix,2);

                % add class coefficient
                u_class=u_nodes;
                if self.sym,u_class=(u_nodes/2)+0.5;end
                fit_matrix_class=fit_matrix.*baseFcnClass(u_class,self.N1,self.N2);

                % reverse calculate control point
                poles=zeros(pole_num,size(nodes,2));
                dim_idx=1:size(nodes,2);dim_idx(2)=[];
                poles(:,dim_idx)=fit_matrix\nodes(:,dim_idx);
                poles(:,2)=fit_matrix_class\nodes(:,2);

                % updata coefs
                weights=ones(pole_num,1);
                poles=[poles.*weights,weights];
                self.coefs=poles;
            end

            pnt_list=self.calPoint(u_nodes);
            RMSE=sqrt(sum((nodes-pnt_list).^2,"all"));
        end
    end

    methods % calculate coord
        function u_list=calCoordinate(self,pnts_init,geom_torl)
            % base on X, Y, Z calculate local coordinate in curve
            %
            if nargin < 3, geom_torl=[];end
            if isempty(geom_torl), geom_torl=sqrt(eps);end

            % find point to start
            u_list=self.findNearest(pnts_init,20);

            % use project function to adjust parameter
            u_list=self.projectPoint(pnts_init,geom_torl,u_list);
        end

        function u_list=projectPoint(self,pnts_init,geom_torl,u_list)
            % adjust U by Jacobian transformation
            % also can project point to curve
            %
            % input:
            % pnt_list(matrix): point_number x dimension matrix
            %
            if nargin < 4
                u_list=[];
                if nargin < 3
                    geom_torl=[];
                end
            end
            if isempty(geom_torl), geom_torl=sqrt(eps);end
            self=self.deriv(1);

            % find point to start
            if isempty(u_list)
                u_list=self.findNearest(pnts_init,20);
            end
            [pnt_num,~]=size(pnts_init);u_list=u_list(:);
            u_min=0+geom_torl;u_max=1-geom_torl;

            % iteration
            iter=0;iter_max=50;
            done=false;
            pnt_idx=1:pnt_num;
            while ~done
                [pnts,dpnts_du]=self.calGradient(u_list(pnt_idx));
                dpnts=pnts_init(pnt_idx,:)-pnts;

                % Jacobian transformation
                RU_RU=sum(dpnts_du.*dpnts_du,2);
                RU_D=sum(dpnts_du.*dpnts,2);
                dus=RU_D./RU_RU;
                dus(isnan(dus) | isinf(dus))=0;

                u_list(pnt_idx)=u_list(pnt_idx)+dus;
                u_list=max(u_list,u_min);u_list=min(u_list,u_max);

                pnt_idx=pnt_idx((abs(RU_D) > geom_torl));

                % pnts_inv=self.calPoint(u_list);
                % scatter3(pnts_inv(:,1),pnts_inv(:,2),pnts_inv(:,3));

                iter=iter+1;
                if isempty(pnt_idx) || iter >= iter_max
                    done=true;
                end
            end
        end

        function u_list=findNearest(self,pnt_list,param)
            % find nearest U in grid
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
            pnt_base_list=self.calPoint(u_base_list);
            u_list=zeros(pnt_num,1);
            for pnt_idx=1:pnt_num
                pnt=pnt_list(pnt_idx,:);
                dis_abs=vecnorm(pnt_base_list-pnt,2,2);
                [~,idx]=min(dis_abs,[],"all");
                u_list(pnt_idx)=u_base_list(idx(1));
            end
        end
    end
end

%% function

function coef_list=convertBezierToPoly(pole_list)
% convert bezier pole to standard polynomial
%
deg=size(pole_list,1)-1;
fac_list=[1,cumprod(1:deg)]';

ir_mat=zeros(deg+1);
sign_mat=zeros(deg+1);
i=0:deg;
ni_list=fac_list(deg+1)./fac_list(deg-i+1)./fac_list(i+1);
for i=0:deg
    r=0:i;
    ir_mat(i+1,1:(i+1))=fac_list(i+1)./fac_list(i-r+1)./fac_list(r+1);
    sign_mat(i+1,1:(i+1))=(-1).^(i-r);
end

coef_list=(ni_list.*ir_mat.*sign_mat)*pole_list;
end

function pole_list=convertPolyToBezier(coef_list)
% convert standard polynomial to bezier pole
%
deg=size(coef_list,1)-1;
fac_list=[1,cumprod(1:deg)]';

miij_mat=zeros(deg+1);
i=0:deg;
mi_list=fac_list(deg+1)./fac_list(deg-i+1)./fac_list(i+1);
for i=0:deg
    j=0:i;
    miij_mat(i+1,1:(i+1))=fac_list(deg-j+1)./fac_list((deg-j)-(i-j)+1)./fac_list((i-j)+1);
end

pole_list=(miij_mat./mi_list)*coef_list;
end

function [pnts,wts]=calCurve(~,number,coefs,u_knots,u_order,us)
% calculate spline point and weight
%
us=us(:);

if number == 2
    pnts=coefs(1,:).*(1-us)+coefs(2,:).*us;
else
    min_u=u_knots(1);max_u=u_knots(end);
    us_k=us*(max_u-min_u)+min_u;

    [u_N_list,u_idx_srt,~]=baseFcnN(us_k,u_order,u_knots);

    % evaluate along the u direction
    pnts=zeros(length(us),size(coefs,2));
    for deg_idx=1:u_order+1
        pnts=pnts+u_N_list(:,deg_idx).*coefs(u_idx_srt+(deg_idx-1),:);
    end
end

wts=pnts(:,end);
pnts=pnts(:,1:end-1);
if nargout < 2
    pnts=pnts./wts;
end
end

function [pnts,dpnts_du]=differFcn1Robust(fcn,us,step)
% calculate function gradient by robust differ
%
pnts=fcn(us);
pnts_uf=fcn(min(us+step,1));
pnts_ub=fcn(max(us-step,0));
bools_f=(us+step) >= 1;
bools_b=(us-step) <= 0;
bools_c=~any([bools_f,bools_b],2);
dpnts_du=zeros(size(pnts)); % allocate memory
dpnts_du(bools_c,:)=(pnts_uf(bools_c,:)-pnts_ub(bools_c,:))/2/step;
dpnts_du(bools_f,:)=(pnts(bools_f,:)-pnts_ub(bools_f,:))/step;
dpnts_du(bools_b,:)=(pnts_uf(bools_b,:)-pnts(bools_b,:))/step;
dpnts_du=real(dpnts_du);
end
