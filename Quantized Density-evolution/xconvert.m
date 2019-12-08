function [y,ofl,ufl] = xconvert(bins, func, bin_convert, result_bins, coeff)

% idea: map blockwise from func, defined over bins, into result_bins
%
% bin_convert contains the converted bin boundaries -- 
% bin_convert(1,:) contains mapping of lower bin limit,
% bin_convert(2,:) contains mapping of upper bin limit.

t_result_bins = result_bins(1):result_bins(2):(result_bins(1) + result_bins(2)*(result_bins(3)-1));
result_func = zeros(1,result_bins(3));

% keep track of overflow/underflow
ofl = 0;
ufl = 0;

for c = 1:bins(3)

	% remember: t_result_bins gives the bin *centers*

	rounds = round((bin_convert(:,c) - result_bins(1))/result_bins(2)) + 1;

	% first: are the indices in range?

	flag = 0;
	partflag = 0;
	min_result_bin = result_bins(1) - 0.5*result_bins(2);
	max_result_bin = result_bins(1) + result_bins(2)*(result_bins(3)-1) + 0.5*result_bins(2);

	% lower range exceeded by both
	if ((rounds(1) < 1) & (rounds(2) < 1))
		ufl = ufl + func(c)*bins(2);
		flag = 1;
	end

	% higher range exceeded by both
	if ((rounds(1) > result_bins(3)) & (rounds(2) > result_bins(3)))
		ofl = ofl + func(c)*bins(2);
		flag = 1;
	end

	% lower range exceeded by a single index
	if ((flag == 0) & (rounds(1) < 1))
		rounds(1) = 1;
		bin_convert(1,c) = min_result_bin;
		partflag = 1;
	end

	% higher range exceeded by a single index
	if ((flag == 0) & (rounds(2) > result_bins(3)))
		rounds(2) = result_bins(3);
		bin_convert(2,c) = max_result_bin;
		partflag = 1;
	end

	if (flag == 0)

		% is the conversion restricted to a single bin?

		if (rounds(1) == rounds(2))

			% map the entire probability of the bin to a single result bin

			bar = func(c)*bins(2)/result_bins(2);
			result_func(rounds(1)) = result_func(rounds(1)) + bar;

		% is the conversion over two bins only?

		elseif ((rounds(2)-rounds(1)) == 1)

			% find the fractional probabilities associated with each bin
			% the bin boundary is (obviously) halfway between round(2) and round(1)

			bdy = (((rounds(1)-1)+0.5)*result_bins(2) + result_bins(1));
			lowfrac = abs(bin_convert(1,c) - bdy)/result_bins(2);
			highfrac = abs(bin_convert(2,c) - bdy)/result_bins(2);
			bar = [lowfrac highfrac] * func(c);
			bar = bar .* [coeff(rounds(1)) coeff(rounds(2))];
			result_func(rounds(1):rounds(2)) = result_func(rounds(1):rounds(2)) + bar;

		else

			% find the fractional probabilities associated with the end bins
			% then the probabilities associated with the intervening bins

			lowbdy = (((rounds(1)-1)+0.5)*result_bins(2) + result_bins(1));
			highbdy = (((rounds(2)-1)-0.5)*result_bins(2) + result_bins(1));
			lowfrac = abs(bin_convert(1,c) - lowbdy)/result_bins(2);
			highfrac = abs(bin_convert(2,c) - highbdy)/result_bins(2);

			bar = zeros(1,rounds(2)-rounds(1)+1);
			bar(1) = lowfrac*func(c);
			bar(rounds(2)-rounds(1)+1) = highfrac*func(c);
			bar(2:(rounds(2)-rounds(1))) = func(c);
			bar = bar .* coeff(rounds(1):rounds(2));
			result_func(rounds(1):rounds(2)) = result_func(rounds(1):rounds(2)) + bar;

		end % if

	end % if

	if (partflag == 1)

		% part of the probability lies outside of the range
		% calculate how much is accounted for; rest is overflow

		pprob = sum(bar)*result_bins(2);

		if (pprob < func(c)*bins(2))

			if (rounds(1) < 1) % underflow
				ufl = ufl + (func(c)*bins(2) - pprob);
			else
				ofl = ofl + (func(c)*bins(2) - pprob);
			end

		end

	end

end % for

y = result_func;
 
