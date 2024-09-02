function Point=meshMapGrid(line_u0,line_u1,line_0v,line_1v,u_list,v_list)
% generate discrete area by mapping
% using Lagrange polynomial mapping method
% local parameter u and v is equispaced
%
% input:
% line_u0 (matrix): u_num x dimension
% line_u1 (matrix): u_num x dimension
% line_v0 (matrix): v_num x dimension
% line_v1 (matrix): v_num x dimension
%
% output:
% Point (matrix): v_num x u_num x dimension
%
if nargin < 6,v_list=[];if nargin < 5,u_list=[];end;end

geom_torl=100*eps;
line_list={line_u0,line_1v,flipud(line_u1),flipud(line_0v)};
line_list=GeomApp.correctLine(line_list,geom_torl);
if ~any(norm(line_list{4}(end,:)-line_list{1}(1,:)) < geom_torl)
    error('GeomApp.MapGrid: line not connect');
end

line_u0=line_list{1};
line_1v=line_list{2};
line_u1=flipud(line_list{3});
line_0v=flipud(line_list{4});

u_num=size(line_u0,1);v_num=size(line_0v,1);
if u_num ~= size(line_u1,1) || v_num ~= size(line_1v,1)
    error('GeomApp.MapGrid: opposite line of boundary discrete number no equal')
end
dimension=size(line_u0,2);
if isempty(u_list),u_list=linspace(0,1,u_num);end,u_list=u_list(:)';
if isempty(v_list),v_list=linspace(0,1,v_num);end,v_list=v_list(:)';

Point=zeros(v_num,u_num,dimension);

% preproces boundary
Point(1,:,:)=line_u0(:,:);
Point(end,:,:)=line_u1(:,:);
Point(:,1,:)=line_0v(:,:);
Point(:,end,:)=line_1v(:,:);

if u_num > 2 && v_num > 2
    H=diff(u_list);G=diff(v_list);

    % solve prepare
    H_ipj=H(1:end-1)+H(2:end);
    H_ij=H(1:end-1).*H(2:end).*H_ipj/2;
    G_ipj=G(1:end-1)+G(2:end);
    G_ij=G(1:end-1).*G(2:end).*G_ipj/2;

    % D_xx
    d_xx=spdiags([H(1:end-1)./H_ij;-H_ipj./H_ij;H(2:end)./H_ij]',-1:1,u_num-2,u_num-2)';
    I_M=speye(v_num-2,v_num-2);
    D_xx=kron(d_xx,I_M);

    % D_yy
    d_yy=spdiags([G(1:end-1)./G_ij;-G_ipj./G_ij;G(2:end)./G_ij]',-1:1,v_num-2,v_num-2)';
    I_K=speye(u_num-2,u_num-2);
    D_yy=kron(I_K,d_yy);

    % Laplace matrix
    Lap=(D_xx+D_yy);

    Point_B=zeros(v_num-2,u_num-2,dimension);
    Point_B(1,:,:)=Point_B(1,:,:)-Point(1,2:u_num-1,:)*(G(2)/G_ij(1));
    Point_B(end,:,:)=Point_B(end,:,:)-Point(end,2:u_num-1,:)*(G(end-1)/G_ij(end));
    Point_B(:,1,:)=Point_B(:,1,:)-Point(2:v_num-1,1,:)*(H(2)/H_ij(1));
    Point_B(:,end,:)=Point_B(:,end,:)-Point(2:v_num-1,end,:)*(H(end-1)/H_ij(end));
    for dim_idx=1:dimension
        Point_B_page=Point_B(:,:,dim_idx);
        Point(2:v_num-1,2:u_num-1,dim_idx)=reshape(Lap\Point_B_page(:),v_num-2,u_num-2);
    end
end
end