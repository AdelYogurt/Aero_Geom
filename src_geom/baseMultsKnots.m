function [mults,knots]=baseMultsKnots(knotvctr)
% base on knot vector to create mults and knots
%
knots=unique(knotvctr);
mults=ones(size(knots));
for k_idx=1:length(knots)
    mults(k_idx)=sum(knotvctr == knots(k_idx));
end
end