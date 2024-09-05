classdef Curve
    % Non-Uniform B-Splines Curve
    % define reference to step standard
    %
    properties % Explicit Attributes
        form='B-NURBS'; % curve_form
        coef_dim=[]; % dimension of coefs
        u_coef_num=[]; % number of coefs
        coefs=[]; % curve coefs matrix, equal to permute(cat(2,poles.*weights,weights),[2,1]);
        u_knotvctr=[]; % knot_vector
        u_order=[]; % degree
    end

    properties % Derivate Attributes
        % step standard properties
        knot_spec=''; % knot_spec
        u_closed=false; % Periodic/closed_curve (boolean)
        intersected=false; % self_intersect (boolean)

        deriv_crv=Curve.empty(0);
    end

    methods % define curve
        function self=Curve(poles,u_degree,u_mults,u_knots,weights)
            % generate Non-Uniform Rational B-Splines Curve
            %
            % calling:
            % crv=Curve(poles,degree,mults,knots,weights)
            % crv=Curve(coefs,knots)
            %
            % input:
            % poles (matrix): control point, pole_num x dimension matrix
            % degree (matrix): optional input
            % mults (matrix): optional input
            % knots (matrix): optional input
            % weights (matrix): optional input
            %
            % output:
            % Curve
            %
            % notice:
            % degree default is pole_num-1, which will be Bezier curve
            %
            if nargin < 5
                weights=[];
                if nargin < 4
                    u_knots=[];
                    if nargin < 3
                        u_mults=[];
                        if nargin < 2
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

                if u_pole_num < (u_order+1)
                    error('Volume: pole_num less than order+1');
                end
            else
                % u default value
                if isempty(u_mults) && isempty(u_knots)
                    u_mults=[u_degree+1,ones(1,u_pole_num-u_degree-1),u_degree+1];
                    u_knots=linspace(0,1,u_pole_num-u_degree+1);
                elseif ~isempty(u_mults) && isempty(u_knots)
                    u_knots=linspace(0,1,length(u_mults));
                elseif isempty(u_mults) && ~isempty(u_knots)
                    error('Curve: need mults input');
                end
                u_knotvrtc=baseKnotVctr(u_mults,u_knots);
                u_knotvrtc=sort(u_knotvrtc);

                if isempty(weights), weights=ones(1,u_pole_num);end
                if ~isempty(weights), weights=weights(:);end

                if u_pole_num < (u_degree+1)
                    error('Curve: pole_num less than degree+1');
                end

                if length(u_knotvrtc) ~= u_pole_num+u_degree+1
                    error('Curve: knot_num is not equal to pole_num+degree+1');
                end

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

        function self=deriv(self,deriv_time)
            % generate derivation curve
            %

            % recursion to derivate data
            if deriv_time > 0 && self.u_order > 1
                if isempty(self.deriv_crv)
                    [coefs_new,~,knots_new]=BSpline.deriv(self.coefs,self.u_order,self.u_knotvctr);
                    self.deriv_crv=Curve(coefs_new,knots_new);
                end

                % donot output data to avoid matlab repeat deriv_crv
                self.deriv_crv=self.deriv_crv.deriv(deriv_time-1);
            end
        end

        function [poles,weights]=getPoles(self)
            % convert coefs to poles and weights
            %
            poles=self.coefs(:,1:end-1);
            weights=self.coefs(:,end);
            weights(weights == 0)=1;
            poles=poles./weights;
        end

        function [pnts,wts]=calPoint(self,us)
            % calculate point on curve
            %
            us=us(:);

            [u_N_list,u_idx_srt,~]=baseFcnN(us,self.u_order,self.u_knotvctr);

            % evaluate along the u direction
            pnts=zeros(length(us),size(self.coefs,2));
            for deg_idx=1:self.u_order+1
                pnts=pnts+u_N_list(:,deg_idx).*self.coefs(u_idx_srt+(deg_idx-1),:);
            end

            wts=pnts(:,end);
            pnts=pnts(:,1:end-1);
            if nargout < 2
                pnts=pnts./wts;
            end
        end

        function varargout=calGradient(self,us)
            % calculate gradient of curve
            %
            deriv_time=nargout-1;
            self=self.deriv(deriv_time);

            % recursion to calculate gradient
            if deriv_time > 0
                % calculate initial curve point and weight
                [c0,w0]=self.calPoint(us);
                p0=c0./w0;

                % calculate first derivative
                [c1,w1]=self.deriv_crv.calPoint(us);
                p1=(c1-p0.*w1)./w0;

                varargout={p0,p1};
                if nargout > 2
                    % calculate second derivative
                    [c2,w2]=self.deriv_crv.deriv_crv.calPoint(us);
                    p2=(c2-2*p1.*w1)./w0-p0.*w2;
                    varargout=[varargout,{p2}];

                    if nargout > 3
                        % calculate thrid derivative
                        [c3,w3]=self.deriv_crv.deriv_crv.deriv_crv.calPoint(us);
                        p3=(c3-3*p2.*w1-3*p1.*w2)./w0-p0.*w3;
                        varargout=[varargout,{p3}];

                        if nargout > 4
                            % calculate fourth derivative
                            [c4,w4]=self.deriv_crv.deriv_crv.deriv_crv.deriv_crv.calPoint(us);
                            p4=(c4-4*p3.*w1-6*p2*w2-4*p1.*w3)./w0-p0.*w4;
                            varargout=[varargout,{p4}];
                        end
                    end
                end
            end
        end

        function varargout=displayGeom(self,axe_hdl,crv_option,param)
            % display curve on axes handle
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
            if isempty(crv_option),crv_option=struct();end
            if isempty(param),param=51;end
            if length(param) == 1,param=linspace(self.u_knotvctr(1),self.u_knotvctr(end),param);end
            pnts=self.calPoint(param);

            % draw points on axe_hdl
            if self.coef_dim-1 == 2
                ln_hdl=line(axe_hdl,pnts(:,1),pnts(:,2),crv_option);
            elseif self.coef_dim-1 == 3
                ln_hdl=line(axe_hdl,pnts(:,1),pnts(:,2),pnts(:,3),crv_option);
                zlabel('\itZ');
            end
            xlabel('\itX');
            ylabel('\itY');

            if nargout > 0, varargout={ln_hdl};end
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
                pole_option=struct('Marker','s','LineStyle','--','Color','r');
            end

            % draw poles on axe_hdl
            poles=self.getPoles();
            if self.coef_dim-1 == 2
                ln_hdl=line(axe_hdl,poles(:,1),poles(:,2),pole_option);
            elseif self.coef_dim-1 == 3
                ln_hdl=line(axe_hdl,poles(:,1),poles(:,2),poles(:,3),pole_option);
                zlabel('\itZ');
            end
            xlabel('\itX');
            ylabel('\itY');

            if nargout > 0, varargout={ln_hdl};end
        end

        function varargout=displayDirect(self,axe_hdl)
            % draw curve on figure handle
            %
            if nargin < 2
                axe_hdl=[];
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % calculate coord on curve
            pnt_mesh=self.coefs(1:2,1:end-1);
            wt_mesh=self.coefs(1:2,end);
            wt_mesh(wt_mesh == 0)=1;
            pnt_mesh=pnt_mesh./wt_mesh;
            origin=pnt_mesh(1,:);
            coord=pnt_mesh(2,:)-origin;
            coord=(coord./vecnorm(coord,2,2))';

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
                zlabel('\itZ');
            end
            xlabel('\itX');
            ylabel('\itY');

            if nargout > 0, varargout={quv_hdl};end
        end
    end

    methods % control curve
        function self=reverseU(self)
            % revese direction of curve
            %
            self.coefs=self.coefs(end:-1:1,:);
            min_k=min(self.u_knotvctr);max_k=max(self.u_knotvctr);dk=max_k-min_k;
            self.u_knotvctr=dk-(self.u_knotvctr(end:-1:1)-min_k)+min_k;
        end

        function self=addOrder(self,times)
            % increase curve degree
            %
            if times <= 0, return;end

            % modify order
            [self.coefs,self.u_order,self.u_knotvctr]=BSpline.addOrder(self.coefs,self.u_order,self.u_knotvctr,times);
            self.u_coef_num=size(self.coefs,2);
        end

        function self=insertKnot(self,knot_ins)
            % insert knot to curve
            %
            if isempty(knot_ins), return;end

            % modify u_list, mults, knots and Ctrls
            [self.coefs,self.u_order,self.u_knotvctr]=BSpline.insertKnot(self.coefs,self.u_order,self.u_knotvctr,knot_ins);
            self.u_coef_num=size(self.coefs,1);
        end

        function self=insertDimension(self,dim_add,bias)
            % add additional geometry dimension
            %
            if nargin < 3 || isempty(bias), bias=0;end
            self.coefs=[self.coefs(:,1:dim_add),...
                ones(self.u_coef_num,1)*bias.*self.coefs(:,end),...
                self.coefs(:,(dim_add+1):end)];
            self.coef_dim=self.coef_dim+1;
        end

        function [crv_1,crv_2]=splitCurve(self,knot_b)
            % split curve at ub
            %
            if knot_b <= min(self.u_knotvctr) || knot_b >= max(self.u_knotvctr)
                error('Curve.splitCurve: knot_b out of boundary of u_list');
            end

            % insert ub to modify poles
            rep_tim=sum(self.u_knotvctr == knot_b);
            self=self.insertKnot(repmat(knot_b,1,self.u_order-rep_tim));

            % locate u place
            knot_num_b=find(self.u_knotvctr == knot_b,1,'last');
            pole_num_b=knot_num_b-self.u_order;

            % generate new curve
            coefs_1=self.coefs(1:pole_num_b,:);
            knots_1=[self.u_knotvctr(1:knot_num_b),knot_b];
            crv_1=Curve(coefs_1,knots_1);

            coefs_2=self.coefs(pole_num_b:end,:);
            knots_2=[knot_b,self.u_knotvctr(knot_num_b-self.u_order+1:end)];
            crv_2=Curve(coefs_2,knots_2);
        end

        function self=translate(self,tran_vctr)
            % translate curve
            %
            tran_vctr=reshape(tran_vctr,1,[]);
            self.coefs(:,1:end-1)=self.coefs(:,1:end-1)+self.coefs(:,end).*tran_vctr(1:self.coef_dim-1);
        end

        function self=rotate(self,rot_mat,rot_cntr)
            % rotate curve
            %
            if nargin < 3 || isempty(rot_cntr), rot_cntr=zeros(1,self.coef_dim-1);end
            self=self.translate(-rot_cntr);
            self.coefs(:,1:end-1)=self.coefs(:,1:end-1)*rot_mat';
            self=self.translate(rot_cntr);
        end
    end

    methods % calculate coord
        function u_list=calCoordinate(self,pnts_init,geom_tol)
            % base on X, Y, Z calculate local coordinate in curve
            %
            if nargin < 3, geom_tol=[];end
            if isempty(geom_tol), geom_tol=sqrt(eps);end

            % find point to start
            u_list=self.findNearest(pnts_init,20);

            % use project function to adjust parameter
            u_list=self.projectPoint(pnts_init,geom_tol,u_list);
        end

        function u_list=projectPoint(self,pnts_init,geom_tol,u_list)
            % adjust U by Jacobian transformation
            % also can project point to curve
            %
            % input:
            % pnt_list(matrix): point_number x dimension matrix
            %
            if nargin < 4
                u_list=[];
                if nargin < 3
                    geom_tol=[];
                end
            end
            if isempty(geom_tol), geom_tol=sqrt(eps);end
            self=self.deriv(1);

            % find point to start
            if isempty(u_list)
                u_list=self.findNearest(pnts_init,20);
            end
            [pnt_num,~]=size(pnts_init);u_list=u_list(:);
            u_min=self.u_knotvctr(1);u_max=self.u_knotvctr(end);

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

                pnt_idx=pnt_idx((abs(RU_D) > geom_tol));

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
            u_base_list=base_list*(self.u_knotvctr(end)-self.u_knotvctr(1))+self.u_knotvctr(1);
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
