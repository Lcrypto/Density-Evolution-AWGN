function y = new_xchk(ext, f_ext, num, mapping)

% for some reason, "wrap" = normal log-likelihood ratio
% and "flog" = log(tanh(x/2)), where x is the llr

% excess = [ofl_pos ufl_pos ofl_neg ufl_neg]
[f_log_pos,f_log_neg,excess] = new_wrap2flog(ext, f_ext, mapping);

p_zero = f_ext(round((ext(3)-1)/2 + 1))*ext(2);

% n_log contains the new bins for the convolved function

n_log = [mapping(1)*num mapping(2) ((mapping(3)+1)*num - num + 1)];
[f_n_log_pos,f_n_log_neg,p_result_zero] = new_fft_convolve_chk(f_log_pos, f_log_neg, num, n_log, excess, p_zero);

[f_n_ext,ofl_pos,ofl_neg] = new_flog2wrap(n_log, f_n_log_pos, f_n_log_neg, ext, p_result_zero);

% note: the overflows should be placed at the value corresponding to
% 2*atanh(exp(-mapping(2)/2)), because this is the uppermost value that
% can be represented using the quantization given by mapping

f_n_ext = new_chk_overflow(f_n_ext,ofl_pos,ofl_neg,ext,mapping); 

y = f_n_ext;
