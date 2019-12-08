function y = new_chk_overflow(f_n_ext,ofl_pos,ofl_neg,ext,mapping)

m_pos_index = round((2*atanh(exp(-mapping(2)/2)) - ext(1))/ext(2)) + 1;
m_neg_index = round((-2*atanh(exp(-mapping(2)/2)) - ext(1))/ext(2)) + 1;

if (m_pos_index < ext(3))
	f_n_ext(m_pos_index) = f_n_ext(m_pos_index) + ofl_pos/ext(2);
else
	f_n_ext(ext(3)) = f_n_ext(ext(3)) + ofl_pos/ext(2);
end

if (m_neg_index > 1)
	f_n_ext(m_neg_index) = f_n_ext(m_neg_index) + ofl_neg/ext(2);
else
	f_n_ext(1) = f_n_ext(1) + ofl_neg/ext(2);
	% Thanks to Sang Hyun Lee for spotting this error
end

% actually, the if-then should be unnecessary, by the nature of overflow ... but 
% then we'd have to check whether ofl_pos or ofl_neg was zero

y = f_n_ext;
