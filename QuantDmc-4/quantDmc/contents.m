%QuantDMC - Matlab library for quantizing discrete memoryless channels.
%Version 4
%
%   Quantizing discrete memoryless channels (DMCs) to maximize mutual 
%   information:
%      quantBiDmc       - Gives one optimal quantizer of a binary-input DMC.  
%                         MEX implementation.
%      quantBiDmcMulti  - Gives all optimal quantizers of a binary-input DMC.
%      quantDmcKLmeans  - KL-means quantization of a DMC.
%      quantDmcGreedy   - Greedy combining quantization of a DMC.
%
%   For binary-input channels QUANTBIDMC should be used, if only one optimal
%   quantzer is needed.  Otherwise, use QUANTBIDMCMULTI, which is slower.
%   For non-binary input channels, QUANTDMCKLMEANS runs faster than 
%   QUANTDMCGREEDY.  When the number of quantizer outputs K is small, 
%   QUANTKLMEANS also has better peformance.
%
%   Additional functions:
%      biAwgn2Dmc        - create a DMC by quantizing an AWGN channel
%      randomDMC         - genreates a random DMC
%      joint2MI          - Compute mutual information 
%      jointDistribution - compute joint distribution from P(Y|X) and P(X)
%      channelSort       - Sort the outputs of a binary-input DMC.
%      mexify            - compile QUANTBIDMC MEX function from C.  QuantDMC 
%                          includes MEX files for some platforms, so is 
%                          required only if you are using another platform.
%
%   Examples are in the example/ directory.
%
%   QuantDMC (c) Brian Kurkoski and contributors
%   Distributed under an MIT-like license; see the file LICENSE

