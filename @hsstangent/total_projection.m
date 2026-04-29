function TH = total_projection(TH)
%TOTAL_PROJECTION Project a tangent vector on the total space.

skew = @(X) .5 *(X - X');

% Projections on the Stiefel manifolds
TH.TU = skew(TH.TU);
TH.TV = skew(TH.TV);
TH.TRl = skew(TH.TRl);
TH.TRr = skew(TH.TRr);
TH.TWl = skew(TH.TWl);
TH.TWr = skew(TH.TWr);

if ~TH.leafnode
    TH.TA11 = total_projection(TH.TA11);
    TH.TA22 = total_projection(TH.TA22);
end

end

