r = 2;
hssoption('block-size', 4);
A = hssgallery('rand', 8, r);
TA = hsstangent('randn', A);

% This will not be vertical!!!!
TA = total_projection(random_vertical(TA, r));

[T, V, M, rhs, C] = projection(TA);
K = [M'*M, C' ; C , zeros(size(C,1))];
Check = projection(V);
