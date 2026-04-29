function TH = mrdivide(TH0, alpha)
%MRDIVIDE Right scalar division for tangent vectors.
%
%   TH = TH0 / alpha

if isa(TH0, 'hsstangent') && isscalar(alpha)
    TH = TH0 * (1 / alpha);
else
    error('hsstangent/mrdivide supports only hsstangent / scalar');
end

end
