function TH = minus(TH1, TH2)
%MINUS Subtraction of two tangent vectors.
%
% Only fields starting with T are combined; base-point data is unchanged.

if ~(isa(TH1, 'hsstangent') && isa(TH2, 'hsstangent'))
    error('hsstangent/minus supports only hsstangent - hsstangent');
end

if ~isequal(base_point(TH1), base_point(TH2))
    error('The vectors belong to different tangent spaces');
end

TH = TH1;

if ~TH1.leafnode
	TH.TB12 = TH1.TB12 - TH2.TB12;
	TH.TB21 = TH1.TB21 - TH2.TB21;
	TH.TR = TH1.TR - TH2.TR;
	TH.TW = TH1.TW - TH2.TW;
	TH.PR = TH1.PR - TH2.PR;
	TH.PW = TH1.PW - TH2.PW;
    TH.TA11 = TH1.TA11 - TH2.TA11;
    TH.TA22 = TH1.TA22 - TH2.TA22;
else
	TH.TU = TH1.TU - TH2.TU;
	TH.TV = TH1.TV - TH2.TV;
	TH.TD = TH1.TD - TH2.TD;
	TH.PU = TH1.PU - TH2.PU;
	TH.PV = TH1.PV - TH2.PV;
end

end
