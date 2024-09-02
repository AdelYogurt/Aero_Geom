function [line_list,map_list,reverse_list]=correctLine(line_list,geom_torl)
% correct line point order
% order should be anti clockwise and start from first line
%
line_num=length(line_list);
map_list=1:line_num;
reverse_list=false(line_num,1);

% start from first line
line_idx=1;
while line_idx <= line_num-1
    line_curr=line_list{line_idx};
    line_next=line_list{line_idx+1};
    if norm(line_curr(end,:)-line_next(1,:)) > geom_torl
        % load remain line all vertex
        vertex_list=zeros((line_num-line_idx)*2,size(line_curr,2));
        for remain_idx=line_idx+1:line_num
            line_rema=line_list{remain_idx};
            vertex_list(2*(remain_idx-line_idx)-1,:)=line_rema(1,:);
            vertex_list(2*(remain_idx-line_idx),:)=line_rema(end,:);
        end

        % search next connect point
        dist=vecnorm((line_curr(end,:)-vertex_list),2,2);
        overlap_idx=find(dist < geom_torl,1);
        if ~any(overlap_idx)
            if line_idx == 1 && ~reverse_list(1)
                % check if reverse first line can save
                line_list{1}=flipud(line_list{1});
                reverse_list(1)=true;
                continue;
            else
                error('correctLine: line not connect');
            end
        end
        exchange_idx=ceil(overlap_idx/2)+line_idx;

        % exchange line
        line_temp=line_list{exchange_idx};
        line_list{exchange_idx}=line_next;
        line_list{line_idx+1}=line_temp;
        map_list(line_idx+1)=exchange_idx;

        % reorder point in line order
        if mod(overlap_idx,2) == 0
            line_list{line_idx+1}=flipud(line_list{line_idx+1});
            reverse_list(line_idx+1)=true;
        end
    end

    line_idx=line_idx+1;
end

end
