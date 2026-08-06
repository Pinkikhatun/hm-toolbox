function Y = retraction(H, TH, t)
%RETRACTION Retraction on the HSS quotient manifold.

if nargin < 2
    t = 1;
end

% Project onto the horizontal space
% TH = total_projection(TH);
% TH = projection(TH);

% Base point and HSS rank
% H = base_point(TH);
% r = hssrank(H);

Y = H;

if H.leafnode

    %% Retract U
    dU = H.U * TH.TU + TH.PU;
    Y.U = uf(H.U + t*dU);

    %% Retract V
    dV = H.V * TH.TV + TH.PV;
    Y.V = uf(H.V + t*dV);

    %% Dense block
    Y.D = H.D + t*TH.TD;

else
    % Y.B12 = H.B12 + t*TH.TB12;
    % Y.B21 = H.B21 + t*TH.TB21;
     %% Recursive HSS blocks
    Y.A11 = retraction(H.A11, TH.TA11, t);
    Y.A22 = retraction(H.A11, TH.TA22, t);

    %% Retract transfer matrix R
    R  = [H.Rl; H.Rr];
    dR = R*TH.TR + TH.PR;

    Rnew = uf(R + t*dR);

    m = size(H.Rl,1);
    Y.Rl = Rnew(1:m,:);
    Y.Rr = Rnew(m+1:end,:);

    %% Retract transfer matrix W
    W  = [H.Wl; H.Wr];
    dW = W*TH.TW + TH.PW;

    Wnew = uf(W + t*dW);

    m = size(H.Wl,1);
    Y.Wl = Wnew(1:m,:);
    Y.Wr = Wnew(m+1:end,:);

    %% Core blocks 
    Y.B12 = H.B12 + t*TH.TB12;
    Y.B21 = H.B21 + t*TH.TB21;

end

end

function A = uf(A)
    [L, ~, R] = svd(A, 0); %
    A = L*R';
end