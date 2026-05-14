r = 2;
hssoption('block-size', 4);
A = hssgallery('rand', 8, r);
TA = hsstangent('randn', A);
TA = random_vertical(TA, r);
[T, V, M, rhs] = projection(TA);
