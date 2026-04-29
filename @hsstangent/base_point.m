function H = base_point(TH)
%BASE_POINT 

H = hss;
H.ml = TH.ml; H.mr = TH.mr;
H.nl = TH.nl; H.nr = TH.nr;

H.B12 = TH.B12; H.B21 = TH.B21;
H.U = TH.U; H.V = TH.V;
H.Rl = TH.Rl; H.Rr = TH.Rr;
H.Wl = TH.Wl; H.Wr = TH.Wr;
H.D = TH.D;

H.leafnode = TH.leafnode;
H.topnode = TH.topnode;

if ~TH.leafnode
    H.A11 = base_point(TH.TA11);
    H.A22 = base_point(TH.TA22);
end

end

