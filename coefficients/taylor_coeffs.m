function [c, params] = taylor_coeffs(func, degree, method, params)
assert(strcmp(method, 'taylor'));
assert(isfield(params, 'anchor'));
assert(isfield(params, 'const'));

assert(degree > 1);

anchor = params.anchor;
const = params.const;
c = func(anchor) * (const).^(0:degree) ./ factorial(0:degree);
c = get_monomial_coeffs_at_zero_(c, 1, -anchor);
end

function c2 = get_monomial_coeffs_at_zero_(c, a, b)
degree = length(c)-1;
c = reshape(c, 1, length(c));
c2 = ones(1,degree+1) * sum(c.*(b.^(0:degree)));
for i = 1 : degree
  tmp = 1;
  for j = 0 : i-1
    tmp = tmp .* ((i:degree) - j);
  end
  c2(i+1) = sum(c(i+1:end).*tmp.*(b.^(0:degree-i))*a^i) / factorial(i);
end
end
