function [f_log_pos,f_log_neg,excess] = new_wrap2flog(ext, f_ext, mapping)

% converts a probability from LLR to log(tanh(L/2)) form
% preserving sign information
%
% note, the zero is dropped! ... but not the underflow

t_mapping = mapping(1):mapping(2):(mapping(1) + mapping(2)*(mapping(3)-1));

wrap = [ext(2) ext(2) (ext(3)-1)/2];

foo = zeros(2,wrap(3));
foo(1,:) = (wrap(1) - wrap(2)/2):wrap(2):(wrap(1)+wrap(2)*(wrap(3)-1)- wrap(2)/2);
foo(2,:) = (wrap(1) + wrap(2)/2):wrap(2):(wrap(1)+wrap(2)*(wrap(3)-1) + wrap(2)/2);
foo(1,1) = wrap(1);
foo(2,wrap(3)) = wrap(1)+wrap(2)*(wrap(3)-1);
bin_convert = log(tanh(foo/2));

coeff = 2*exp(t_mapping)./(1-exp(2*t_mapping));

% excess = [ofl_pos ufl_pos ofl_neg ufl_neg]
excess = zeros(1,4);
bar = f_ext(((ext(3)-1)/2 + 2):ext(3));
[f_log_pos,ofl,ufl] = xconvert(wrap, bar, bin_convert, mapping, coeff);
excess(1) = ofl;
excess(2) = ufl;

bar = f_ext(((ext(3)-1)/2):-1:1);
[f_log_neg,ofl,ufl] = xconvert(wrap, bar, bin_convert, mapping, coeff);
excess(3) = ofl;
excess(4) = ufl;

% note that excess(3) should equal (very close to) zero under the symmetry condition
% ..... NO!  It is still nonzero, because of the granularity of the mapping

