function crv=curveGlue(crv_list,geom_torl)
% connect BSpline curve into one curve
%
if nargin < 3,geom_torl=[];end

if isempty(geom_torl),geom_torl=100*eps;end

% search max order
order_tar=0;
crv_idx=1;
while crv_idx <= length(crv_list)
    crv_curr=crv_list(crv_idx);
    if isempty(crv_curr)
        crv_list(crv_idx)=[];
    else
        if order_tar < crv_curr.u_order
            order_tar=crv_curr.u_order;
        end
    end
    crv_idx=crv_idx+1;
end
crv_num=length(crv_list);

% load all poles to correct line order
line_list=cell(1,crv_num);
for crv_idx=1:crv_num
    line_list{crv_idx}=crv_list(crv_idx).coefs;
end
[~,map_list,order_list]=correctLine(line_list,geom_torl);
crv_list=crv_list(map_list);
for crv_idx=1:crv_num
    if order_list(crv_idx),crv_list(crv_idx)=crv_list(crv_idx).reverseU();end
end

% increase order of edg and load data
coefs=[];
u_knotvctr=[];
len_total=0;
for crv_idx=1:crv_num
    % load edge
    crv=crv_list(crv_idx);

    crv=crv.addOrder(order_tar-crv.u_order);
    len_edg=sum(vecnorm(diff(crv.coefs,1,1),2,2));

    coefs=[coefs(1:end-1,:);crv.coefs];
    u_knotvctr=[u_knotvctr(1:end-1-order_tar),crv.u_knotvctr(2:end)*len_edg+len_total];
    len_total=len_total+len_edg;
end
u_knotvctr=[u_knotvctr(1),u_knotvctr];

u_knotvctr=u_knotvctr/len_total;
crv=Curve(coefs,u_knotvctr);
end
