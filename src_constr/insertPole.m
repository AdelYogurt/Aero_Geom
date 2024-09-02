function poles=insertPole(poles,idx)
% insert a new param into poles
%
size_poles=size(poles);
size_poles=[size_poles(1:idx),1,size_poles(idx+1:end)];
poles=reshape(poles,size_poles);
end