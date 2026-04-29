function [TH] = create_randvec(H,TH, r)
%CREATE_RANDVEC Summary of this function goes here
%   Detailed explanation goes here

if ~H.leafnode
%HSSTANGENT Construct an instance of this class
            %   Detailed explanation goes here
            TH.TB12=randn(r);
            TH.TB21=randn(r);
end        
        % Factorization of the upper triangular block as U12 * V12'
        if ~isempty(H.U)
        %collinear component of U and V
            TH.TU = randn(r);
            TH.TU = TH.TU-TH.TU';
            TH.TV = randn(r);
            TH.TV = TH.TV-TH.TV';
        %%Orthogonal component with respect to U and V
            TH.PU  = randn(size(H, 1)-r,r); 
            TH.PV  = randn(size(H, 2)-r,r);
        end    
        % Factorization of the lower triangular block as U21 * V21'
        if ~isempty(H.Rl)
           TH.TRl = randn(r);
           TH.TRl = TH.TRl-TH.TRl';
           TH.PRl = randn(H.ml-r, r);
           TH.TRr = randn(r);
           TH.TRr = TH.TRr-TH.TRr';
           TH.PRr = randn(H.mr-r, r);
           TH.TWl = randn(r);
           TH.TWl = TH.TWl-TH.TWl';
           TH.PWl = randn(H.nl-r, r);
           TH.TWr = randn(r);
           TH.TWr = TH.TWr-TH.TWr';
           TH.PWr =randn(H.nr-r, r);
        end
        % Dense version of the matrix, if the size is smaller than the
        % minimum allowed block size.
        if ~isempty(H.D)
           TH.TD = randn(size(H.D));
        end
        
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
              TH.TA11= hsstangent();
              TH.TA22= hsstangent();
              TH.TA11 = create_randvec(H.A11,TH.TA11, r);
              TH.TA22 = create_randvec(H.A22,TH.TA22, r);
          end
        

end