# Covariance-Evolution

Covariance Evolution Implementation for BEC-channel with Octopus protograph. After paper published, i'll disclosure full version of source code.


Example of waterfall prediction for rate and length adaptive (25/32<Rate<25/74, 10000<k< 250) MET LDPC code (Richardson protograph) under AWGN-channel (use same idea as Yury Polyanskiy Bound, shift from threshold, but not full analytical proof yet, for BSC exist EXIT-like chart with some stronger proof). In Dr. Richardson file tabulated a,b,VN,CN,Circulant and calculate WER(SNR) = Q(a*(10*log10(SNR)-b)). Using linear interpolation get values for circulant between.

