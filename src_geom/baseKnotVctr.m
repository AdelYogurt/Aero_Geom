function knotvctr=baseKnotVctr(mults,knots)
% base on mults and knots to create knot vector
%
if length(mults) ~= length(knots)
    error('baseKnotVctr: length of mults and knots is not equal')
end
knotvctr=zeros(1,sum(mults));
start_idx=1;end_idx=mults(1);
for n_idx=1:length(knots)
    knotvctr(start_idx:end_idx)=knots(n_idx);
    start_idx=start_idx+mults(n_idx);
    end_idx=end_idx+mults(n_idx);
end
end