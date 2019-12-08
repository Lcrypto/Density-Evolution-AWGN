function [f_wrap,ofl_pos,ofl_neg] = new_flog2wrap(n_log, f_n_log_pos, f_n_log_neg, ext, p_n_zero)

% converts a log(tanh(L/2)) form random variable to LLR
% recall that sign information is preserved

wrap = [ext(2) ext(2) (ext(3)-1)/2];
t_wrap = wrap(1):wrap(2):(wrap(1) + wrap(2)*(wrap(3)-1));

foo = zeros(2,n_log(3));
foo(1,:) = (n_log(1) - n_log(2)/2):n_log(2):(n_log(1)+n_log(2)*(n_log(3)-1) - n_log(2)/2);
foo(2,:) = (n_log(1) + n_log(2)/2):n_log(2):(n_log(1)+n_log(2)*(n_log(3)-1) + n_log(2)/2);
foo(1,1) = n_log(1);
foo(2,n_log(3)) = n_log(1) + n_log(2)*(n_log(3)-1);
bin_convert = 2*atanh(exp(foo));

bar = tanh(t_wrap/2);
coeff = (1./bar).*(1-bar.^2)*0.5;

[f_wrap_pos,ofl_pos,ufl_pos] = xconvert(n_log, f_n_log_pos, bin_convert, wrap, coeff);
[f_wrap_neg,ofl_neg,ufl_neg] = xconvert(n_log, f_n_log_neg, bin_convert, wrap, coeff);

result = zeros(1,ext(3));

result(1:((ext(3)-1)/2)) = f_wrap_neg(((ext(3)-1)/2):-1:1);
result(((ext(3)-1)/2 + 2):ext(3)) = f_wrap_pos;
result((ext(3)-1)/2 + 1) = (p_n_zero + ufl_pos + ufl_neg)/ext(2);

f_wrap = result;


