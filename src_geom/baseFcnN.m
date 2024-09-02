function [N_list,idx_srt,idx_end]=baseFcnN(us,k,u_list)
% Basis function for B-Spline
%
% input:
% u: parametric points
% k: spline degree
% u_list: knot sequence
%
% output:
% N_list(matrix): N-Basis functions vector, numel(u)*(k+1)
% idx_srt(vector): calculate start index
% idx_end(vector): calculate end index
%
if any(us < u_list(1)) || any(us > u_list(end))
    error('baseFcnN: u out of boundary of u_list')
end
if any(isnan(us))
    error('baseFcnN: u have nan')
end
us=us(:);

% modify
bool=us <= u_list;
idx_end=max(length(u_list)-sum(bool,2),k+1);
% idx_end=zeros(size(us));
% for j=1:numel(us)
%     if (us(j)==u_list(k+1)),idx_end(j)=k+1; continue,end
%     idx_end(j)=find(us(j) <= u_list,1,'first')-1;
% end
idx_srt=idx_end-k;

u_num=numel(us);
N_list=zeros(u_num,k+1);
% N=zeros(1,k+1); % modify
for u_idx=1:u_num
    i=idx_end(u_idx); %% findspan uses 0-based numbering
    u=us(u_idx);
    left=zeros(k+1,1);
    right=zeros(k+1,1);

    N_list(u_idx,1)=1;
    for j=1:k
        left(j+1)=u-u_list(i+1-j);
        right(j+1)=u_list(i+j)-u;
        saved=0;
        for r=0:j-1
            temp=N_list(u_idx,r+1)/(right(r+2)+left(j-r+1));
            N_list(u_idx,r+1)=saved+right(r+2)*temp;
            saved=left(j-r+1)*temp;
        end
        N_list(u_idx,j+1)=saved;
    end

%     % modify
%     N(1)=1;
%     for j=1:k
%         left(j+1)=u-u_list(i+1-j);
%         right(j+1)=u_list(i+j)-u;
%         saved=0;
%         for r=0:j-1
%             temp=N(r+1)/(right(r+2)+left(j-r+1));
%             N(r+1)=saved+right(r+2)*temp;
%             saved=left(j-r+1)*temp;
%         end
%         N(j+1)=saved;
%     end
%     N_list(u_idx,:)=N;
end
end
