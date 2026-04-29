function TH = plus(TH1, TH2)
%PLUS Addition of two tangent vectors.
%
% Only fields starting with T are combined; base-point data is unchanged.

if ~(isa(TH1, 'hsstangent') && isa(TH2, 'hsstangent'))
    error('hsstangent/plus supports only hsstangent + hsstangent');
end

if ~isequal(base_point(TH1), base_point(TH2))
    error('The vectors belong to different tangent spaces');
end

TH = TH1;

TH.TB12 = TH1.TB12 + TH2.TB12;
TH.TB21 = TH1.TB21 + TH2.TB21;

TH.TU = TH1.TU + TH2.TU;
TH.TV = TH1.TV + TH2.TV;

TH.TRl = TH1.TRl + TH2.TRl;
TH.TRr = TH1.TRr + TH2.TRr;
TH.TWl = TH1.TWl + TH2.TWl;
TH.TWr = TH1.TWr + TH2.TWr;

TH.TD = TH1.TD + TH2.TD;

if ~TH1.leafnode
    TH.TA11 = TH1.TA11 + TH2.TA11;
    TH.TA22 = TH1.TA22 + TH2.TA22;
end

end
