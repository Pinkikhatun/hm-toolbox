function TH = create_randvec_from_cell(H, C, TH)
% top and bottom blocks of the matrix.
TH.B12 = H.B12;
TH.B21 = H.B21;


%Part related HSS matrix
% Factorization of the upper triangular block as U12 * V12'
TH.U = H.U;
TH.V = H.V;

% Factorization of the lower triangular block as U21 * V21'
TH.Rl = H.Rl;
TH.Rr =H.Rr;
TH.Wl = H.Wl;
TH.Wr = H.Wr;

% Size of the matrix
TH.ml =H.ml;
TH.nl =H.nl;
TH.mr =H.mr;
TH.nr =H.nr;

% Dense version of the matrix, if the size is smaller than the
% minimum allowed block size.
TH.D = H.D;

TH.topnode=H.topnode;
TH.leafnode= H.leafnode;

if ~H.leafnode
    TH.TB12 = C{1};
    TH.TB21 = C{2};
    if ~ H.topnode
        TH.TRl = C{3};
        TH.PRl = C{4};
        TH.TRr = C{5};
        TH.PRr = C{6};
        TH.TWl = C{7};
        TH.PWl = C{8};
        TH.TWr = C{9};
        TH.PWr = C{10};
    end
    TH.TA11= hsstangent();
    TH.TA11= create_randvec_from_cell(H.A11, C{end-1}, TH.TA11);
    TH.TA22= hsstangent();
    TH.TA22= create_randvec_from_cell(H.A22, C{end}, TH.TA22);
else
    TH.TU = C{1};
    TH.PU = C{2};
    TH.TV = C{3};
    TH.PV = C{4};
    TH.TD = C{5};
end


end