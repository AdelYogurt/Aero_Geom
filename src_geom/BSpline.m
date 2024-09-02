classdef BSpline
    methods(Static) % BSpline function
        function [coefs_new,u_order_new,u_knotvrtc_new]=addOrder(coefs,u_order,u_knotvrtc,times)
            % increase curve degree
            %
            if times <= 0
                coefs_new=coefs;
                u_order_new=u_order;
                u_knotvrtc_new=u_knotvrtc;
                return;
            end

            [u_mults,u_knots]=baseMultsKnots(u_knotvrtc);

            while times > 0
                % add repeat node
                u_order_new=u_order+1;
                u_mults_new=u_mults+1;
                u_knotvrtc_new=baseKnotVctr(u_mults_new,u_knots);

                % calculate new ctrl
                matrix=zeros(size(coefs,1)+length(u_knots)-1,size(coefs,1));
                for j=1:size(coefs,1)+length(u_knots)-1
                    for i=1:size(coefs,1)
                        matrix(j,i)=baseFcnL(u_knotvrtc,u_knotvrtc_new,i,j,u_order);
                    end
                end
                matrix=matrix/u_order_new;
                coefs_new=matrix*coefs;

                % sort old curve data
                coefs=coefs_new;
                u_order=u_order_new;
                u_mults=u_mults_new;
                u_knotvrtc=u_knotvrtc_new;

                times=times-1;
            end
        end

        function [coefs_new,order_new,knotvrtc_new]=insertKnot(coefs,order,knotvrtc,knot_ins)
            % insert knot to Spline
            %

            % allocate memory
            [coef_num,coef_dim]=size(coefs);
            knot_ins=sort(knot_ins);
            nu=length(knot_ins);
            u_num=length(knotvrtc);
            coef_num_new=coef_num+nu;
            coefs_new=zeros(coef_num_new,coef_dim);
            knotvrtc_new=zeros(1,u_num+nu);

            n=coef_num-1;
            r=nu-1;
            m=n+order+1;

            % locate add knots position
            if (knot_ins(1) == knotvrtc(coef_num+1)),a=coef_num-1;
            else,a=find(knot_ins(1) >= knotvrtc,1,'last')-1;end
            if (knot_ins(r+1) == knotvrtc(coef_num+1)),b=coef_num-1;
            else,b=find(knot_ins(r+1) >= knotvrtc,1,'last')-1;end
            b=b+1;

            coefs_new(1:a-order+1,:)=coefs(1:a-order+1,:);
            coefs_new(b+nu:coef_num+nu,:)=coefs(b:coef_num,:);

            order_new=order;

            knotvrtc_new(1:a+1)=knotvrtc(1:a+1);
            knotvrtc_new(b+order+nu+1:m+nu+1)=knotvrtc(b+order+1:m+1);

            ii=b+order-1;
            ss=ii+nu;

            % calculate new ctrl point
            for jj=r:-1:0
                idx=(a+1):ii;
                idx=idx(knot_ins(jj+1) <= knotvrtc(idx+1));
                coefs_new(idx+ss-ii-order,:)=coefs(idx-order,:);
                knotvrtc_new(idx+ss-ii+1)=knotvrtc(idx+1);
                ii=ii-length(idx);
                ss=ss-length(idx);
                coefs_new(ss-order,:)=coefs_new(ss-order+1,:);
                for l=1:order
                    idx=ss-order+l;
                    alfa=knotvrtc_new(ss+l+1)-knot_ins(jj+1);
                    if abs(alfa) == 0
                        coefs_new(idx,:)=coefs_new(idx+1,:);
                    else
                        alfa=alfa/(knotvrtc_new(ss+l+1)-knotvrtc(ii-order+l+1));
                        tmp=(1-alfa)*coefs_new(idx+1,:);
                        coefs_new(idx,:)=alfa*coefs_new(idx,:)+tmp;
                    end
                end
                knotvrtc_new(ss+1)=knot_ins(jj+1);
                ss=ss-1;
            end
        end
    
        function [coefs_new,order_new,knotvrtc_new]=deriv(coefs,order,knotvrtc)
            % generate BSpline derivative curve
            %
            [coef_num,coef_dim]=size(coefs);

            % calculate deriv coefs
            coefs_new=zeros(coef_num-1,coef_dim);
            for ctrl_idx=0:coef_num-2
                coefs_new(ctrl_idx+1,:)=order/(knotvrtc(ctrl_idx+order+2)-knotvrtc(ctrl_idx+2))...
                    *(coefs(ctrl_idx+2,:)-coefs(ctrl_idx+1,:));
            end

            % decrease order
            order_new=order-1;
            knotvrtc_new=knotvrtc(2:end-1);
        end
    end
end