function TH = mldivide(alpha, TH0)
%MLDIVIDE Left scalar division for tangent vectors.
%
%   TH = alpha \ TH0

if isscalar(alpha) && isa(TH0, 'hsstangent')
    TH = (1 / alpha) * TH0;
else
    error('hsstangent/mldivide supports only scalar \\ hsstangent');
end

end
