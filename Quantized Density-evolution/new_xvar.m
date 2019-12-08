function y = new_xvar(chan, func, num, state, ext)

sc = length(chan);
sf = length(func);

zeropad = zeros(1,sc + sf*num - num);
bar = zeropad;

% direct probabilities
bar(1:sc) = chan*state(2);
zeropad(1:sf) = func*ext(2);

F_zeropad = fft(zeropad);
foo = F_zeropad.^num;
foo = foo .* fft(bar);
IF_zeropad = ifft(foo);

% extract function of the appropriate length from the middle ...
% question: where is the zero?
% clearly, the minimum index is the sum of all input minimum indices

minx = state(1) + num*ext(1);

% let's hope that state(2) and ext(2) are the same

ext_minx_index = round((ext(1) - minx)/ext(2)) + 1;
ext_maxx_index = ext_minx_index + ext(3) - 1;

%thud = IF_zeropad/sum(IF_zeropad);
ufl = abs(sum(IF_zeropad(1:(ext_minx_index-1))));
ofl = abs(sum(IF_zeropad((1+ext_maxx_index):length(IF_zeropad))));

IF_zeropad(ext_minx_index) = IF_zeropad(ext_minx_index) + ufl;
IF_zeropad(ext_maxx_index) = IF_zeropad(ext_maxx_index) + ofl;

y = abs(IF_zeropad(ext_minx_index:ext_maxx_index))/ext(2);


