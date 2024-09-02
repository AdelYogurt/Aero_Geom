function [srf,u_nodes,v_nodes]=interpPointToSurface(nodes,u_degree,v_degree,u_pole_num,v_pole_num,u_nodes,v_nodes)
% generate BSpline curve by defined fitting points
%
% input:
% nodes (matrix): fit point,u_node_num x v_node_num x dimension matrix
% degree (optional):
% mults (optional):
% knots (optional):
% u_pole_num (optional):
% v_pole_num (optional):
% u_nodes (optional):
% v_nodes (optional):
%
% output:
% srf (surface): Spline surface
% u_nodes (matrix): local parameter of interpolated point
% v_nodes (matrix): local parameter of interpolated point
%
% notice:
% Input degree default is pole_num-1, which is Bezier surface
%
if nargin < 7
    v_nodes=[];
    if nargin < 6
        u_nodes=[];
        if nargin < 5
            v_pole_num=[];
            if nargin < 4
                u_pole_num=[];
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

[v_node_num,u_node_num,dimension]=size(nodes);
if isempty(u_pole_num),u_pole_num=u_node_num;end
if isempty(v_pole_num),v_pole_num=v_node_num;end
if u_pole_num > u_node_num || v_pole_num > v_node_num
    error('interpPointToSurface: pole_num more than node_num')
end

% default value of u_nodes
if isempty(u_nodes)
    u_nodes=vecnorm(nodes(2:end,:,:)-nodes(1:end-1,:,:),2,3);
    u_nodes=mean(u_nodes,2);u_nodes=[0;cumsum(u_nodes)];
end
if all(size(u_nodes) ~= 1), error('interpPointToSurface: u_nodes can not be matrix');end
u_nodes=(u_nodes(:)-min(u_nodes))/(max(u_nodes)-min(u_nodes));
u_node_knots=interp1(linspace(0,1,length(u_nodes)),u_nodes,linspace(0,1,u_pole_num));

% default value of v_nodes
if isempty(v_nodes)
    v_nodes=vecnorm(nodes(:,2:end,:)-nodes(:,1:end-1,:),2,3);
    v_nodes=mean(v_nodes,1);v_nodes=[0;cumsum(v_nodes')];
end
if all(size(v_nodes) ~= 1), error('interpPointToSurface: v_nodes can not be matrix');end
v_nodes=(v_nodes(:)-min(v_nodes))/(max(v_nodes)-min(v_nodes));
v_node_knots=interp1(linspace(0,1,length(v_nodes)),v_nodes,linspace(0,1,v_pole_num));

% calculate u_list
u_mults=[u_degree+1,ones(1,u_pole_num-u_degree-1),u_degree+1];
u_knots=linspace(0,1,u_pole_num-u_degree+1);
for j=2:u_pole_num-u_degree
    u_knots(j)=mean(u_node_knots(j:j+u_degree-1));
end
% modify
u_knots=interp1(linspace(0,1,length(u_knots)),u_knots,linspace(0,1,u_pole_num-u_degree+1));
u_list=baseKnotVctr(u_mults,u_knots);

% calculate v_list
v_mults=[v_degree+1,ones(1,v_pole_num-v_degree-1),v_degree+1];
v_knots=linspace(0,1,v_pole_num-v_degree+1);
for j=2:v_pole_num-v_degree
    v_knots(j)=mean(v_node_knots(j:j+v_degree-1));
end
% modify
v_knots=interp1(linspace(0,1,length(v_knots)),v_knots,linspace(0,1,v_pole_num-v_degree+1));
v_list=baseKnotVctr(v_mults,v_knots);

% base on node point list inverse calculate control point list
u_fit_matrix=zeros(u_node_num,u_pole_num);
[N_list,idx_srt,idx_end]=baseFcnN(u_nodes,u_degree,u_list);
for deg_idx=1:u_degree+1
    idx=sub2ind([u_node_num,u_pole_num],(1:(u_node_num))',idx_srt+(deg_idx-1));
    u_fit_matrix(idx)=N_list(:,deg_idx);
end
v_fit_matrix=zeros(v_node_num,v_pole_num);
[N_list,idx_srt,idx_end]=baseFcnN(v_nodes,v_degree,v_list);
for deg_idx=1:v_degree+1
    idx=sub2ind([v_node_num,v_pole_num],(1:(v_node_num))',idx_srt+(deg_idx-1));
    v_fit_matrix(idx)=N_list(:,deg_idx);
end
v_fit_matrix=v_fit_matrix';
poles=pagemrdivide(pagemldivide(u_fit_matrix,nodes),v_fit_matrix);

% adjust bound poles
poles(:,1,:)=pagemldivide(u_fit_matrix,nodes(:,1,:));
poles(:,end,:)=pagemldivide(u_fit_matrix,nodes(:,end,:));
poles(1,:,:)=pagemrdivide(nodes(1,:,:),v_fit_matrix);
poles(end,:,:)=pagemrdivide(nodes(end,:,:),v_fit_matrix);

poles(1,1,:)=nodes(1,1,:);
poles(end,1,:)=nodes(end,1,:);
poles(1,end,:)=nodes(1,end,:);
poles(end,end,:)=nodes(end,end,:);

srf=Surface(poles,u_degree,v_degree,u_mults,v_mults,u_knots,v_knots);
end
