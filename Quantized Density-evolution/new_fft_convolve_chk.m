function [f_n_log_pos,f_n_log_neg,p_result_zero] = new_fft_convolve_chk(f_log_pos, f_log_neg, num, n_log, excess, p_zero)

% perform FFTs over R x GF(2) to obtain probabilities with sign information

% I will *assume* that the range of f_log_pos and f_log_neg are over the range
% min:increment:-increment, with min*increment elements
%
% the input ofl_pos is the probability of the input "flog" message
% being zero ... can directly include this in the sum
%
% furthermore, by the symmetry property, "flog" can only be zero
% if the sign is positive

% find the probabilities of being positive or negative ... used in calculating the
% conditional probability

% these form the extended (by 1) messages, which give room for an extra zero message
q_pos = zeros(1,length(f_log_pos)+1);
q_pos(1:length(f_log_pos)) = f_log_pos;
q_neg = zeros(1,length(f_log_neg)+1); 
q_neg(1:length(f_log_neg)) = f_log_neg;

% sf = size of the extended message
sf = length(f_log_pos)+1;

% find length to nearest higher power of 2 to help with fft
zp2_size = 2^(ceil(log2(sf*num - num + 1)));
zeropad_pos = zeros(1, zp2_size);
zeropad_neg = zeropad_pos;

% preserve *direct* probabilities under convolution
% also note that overflow is included in zeropad_pos
zeropad_pos(1:sf) = q_pos*n_log(2); 
zeropad_neg(1:sf) = q_neg*n_log(2);
zeropad_pos(sf) = excess(1);
zeropad_neg(sf) = excess(3);

% find the probabilities of being positive or negative and finite ... 
% as well as marginal and conditional probabilities
p_pos_fin = sum(zeropad_pos);
p_pos = sum(zeropad_pos) + excess(2);

% note: excess(1), excess(3) are already contained in zeropad

p_neg_fin = sum(zeropad_neg);
p_neg = sum(zeropad_neg) + excess(4);
% fprintf('%f \n ',zeropad_pos/p_pos_fin);
% fprintf('%f \n ',zeropad_neg/p_pos_fin);
% here we take the FFT of the *conditional* density
F_pos_fin = fft(zeropad_pos/p_pos_fin);
F_neg_fin = fft(zeropad_neg/p_neg_fin);

result_pos = zeros(1, zp2_size);
result_neg = zeros(1, zp2_size);
p_result_zero = 0;

% now combine the magnitude and sign components to get the result
% have to do proper handing for -\infty messages (i.e., wrap=0)

foo = zeros(num+1,zp2_size);

% now marginalize with respect to positive or negative events
% (i.e., c even or odd)
%
% note that result_pos and result_neg are joint probabilities
% of the amplitude and being positive or negative

for c = 0:num

	% no underflow
	% any sign = zero -- attached to p_result_zero
	v = (1 - (p_pos/(p_pos+p_zero))^(num-c)) * nchoosek(num,c)*p_neg^c * (p_pos+p_zero)^(num-c);
	p_result_zero = p_result_zero + v;
   temp2=(1-excess(4))^c * (1-excess(2))^(num-c);
   temp3=nchoosek(num,c) * p_neg^c * p_pos^(num-c);
	v = (F_pos_fin.^(num-c)).*(F_neg_fin.^c);
	v = v * temp2;
	v = v * temp3;

	if (mod(c,2)==0)
		% even negatives -- result is positive
		result_pos = result_pos + v;
	else
		% odd negatives -- result is negative
		result_neg = result_neg + v;
	end

	% at least one underflow
	% regardless of sign, attached to p_result_zero

	% if there are c negative messages, then there are at most num-c positive messages

	for d = 0:(num-c)

		% calculate no underflow first ... then complement
		w = (1-excess(2))^c * (1-excess(4))^d;
		v = (1-w)*prod(1:num)/(prod(1:c)*prod(1:d)*prod(1:(num-c-d)));
		v = v * p_neg^c * p_pos^d * p_zero^(num-c-d);

		p_result_zero = p_result_zero + v;

	end

end 
result_pos=ifft(result_pos);
result_neg=ifft(result_neg);
result_pos = abs(result_pos)/n_log(2);
result_neg = abs(result_neg)/n_log(2);

f_n_log_pos = result_pos(1:(sf*num - num + 1));
f_n_log_neg = result_neg(1:(sf*num - num + 1));


