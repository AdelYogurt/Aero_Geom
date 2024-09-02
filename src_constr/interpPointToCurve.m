function [crv,u_nodes]=interpPointToCurve(nodes,u_degree,u_pole_num,u_nodes,derivs)
% generate BSpline curve by defined fitting points
%
% input:
% nodes (matrix): fit point, dimension x node_num matrix
% degree (optional):
% mults (optional):
% knots (optional):
% pole_num (optional):
% u_nodes (optional):
% Derivs (optional): derivative of nodes. Set nan if der is free
%
% output:
% crv (Curve): Spline curve
% u_nodes (matrix): local parameter of interpolated point
%
% notice:
% Input degree default is pole_num-1, which is Bezier curve
%
if nargin < 5
    derivs=[];
    if nargin < 4
        u_nodes=[];
        if nargin < 3
            u_pole_num=[];
            if nargin < 2
                u_degree=[];
            end
        end
    end
end

% derivative number
if isempty(derivs) || all(any(isnan(derivs),2)), u_deriv_num=0;
else, u_deriv_num=sum(all(~isnan(derivs),2));end

[u_node_num,coef_dim]=size(nodes);
if isempty(u_pole_num),u_pole_num=u_node_num+u_deriv_num;end
if isempty(u_degree),u_degree=u_pole_num-1;end
if u_pole_num > u_node_num+u_deriv_num
    error('interpPointToCurve: pole_num more than node_num')
end

% default value of u_nodes vector
if isempty(u_nodes)
    u_nodes=vecnorm(nodes(2:end,:)-nodes(1:end-1,:),2,2);
    u_nodes=[0;cumsum(u_nodes)];
end
u_nodes=(u_nodes(:)-min(u_nodes))/(max(u_nodes)-min(u_nodes));
u_node_knots=u_nodes;

% process derivs of u_nodes vector
if u_node_num ~= 0
    u_node_knots=[u_nodes;u_nodes(all(~isnan(derivs),2))];
    u_node_knots=sort(u_node_knots);
end
u_node_knots=interp1(linspace(0,1,length(u_node_knots)),u_node_knots,linspace(0,1,u_pole_num));

if ~isempty(derivs) && size(derivs,1) ~= u_node_num
    error('interpPointToCurve: Derivs matrix do not equal to node_num');
end

% reference:
% [1] 施法中. 计算机辅助几何设计与非均匀有理B样条[M]. 208-281.
% [2] https://blog.csdn.net/he_nan/article/details/107134271
u_mults=[u_degree+1,ones(1,u_pole_num-u_degree-1),u_degree+1];
u_knots=linspace(0,1,u_pole_num-u_degree+1);
for j=2:u_pole_num-u_degree
    u_knots(j)=mean(u_node_knots(j:j+u_degree-1));
end
u_knots=interp1(linspace(0,1,length(u_knots)),u_knots,linspace(0,1,u_pole_num-u_degree+1)); % modify
u_list=baseKnotVctr(u_mults,u_knots);

% base on node point list inverse calculate control point list
fit_matrix=zeros(u_node_num+u_deriv_num,u_pole_num);
[N_list,idx_srt,idx_end]=baseFcnN(u_nodes,u_degree,u_list);
for deg_idx=1:u_degree+1
    idx=sub2ind([u_node_num+u_deriv_num,u_pole_num],(1:(u_node_num))',idx_srt+(deg_idx-1));
    fit_matrix(idx)=N_list(:,deg_idx);
end

% process derivative of fit matrix
if u_deriv_num ~= 0
    u_list_deriv=u_list(2:end-1);
    U_node_deriv=u_nodes(all(~isnan(derivs),2));
    for deriv_idx=1:u_deriv_num
        u=U_node_deriv(deriv_idx);

        d=u_degree/(u_list(1+u_degree+1)-u_list(1+1));
        fit_matrix(u_node_num+deriv_idx,1)=-d*baseFcnN(u,u_degree-1,u_list_deriv);
        for ctrl_idx=2:u_pole_num-1
            dn=u_degree/(u_list(ctrl_idx+u_degree+1)-u_list(ctrl_idx+1));
            fit_matrix(u_node_num+deriv_idx,ctrl_idx)=...
                d*baseFcnN(u_list_deriv,u,ctrl_idx-1,u_degree-1)-...
                dn*baseFcnN(u_list_deriv,u,ctrl_idx,u_degree-1);
            d=dn;
        end
        fit_matrix(u_node_num+deriv_idx,u_pole_num)=dn*baseFcnN(u_list_deriv,u,u_pole_num-1,u_degree-1);
    end
end

if u_deriv_num == 0, poles=fit_matrix\nodes;
else, poles=fit_matrix\[nodes;derivs(all(~isnan(derivs),2),:)];end

% adjust bound poles
poles(1,:)=nodes(1,:);
poles(end,:)=nodes(end,:);

crv=Curve(poles,u_degree,u_mults,u_knots);
end
