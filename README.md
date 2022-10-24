# Density-Evolution-AWGN
1. main.m at subfolder https://github.com/Lcrypto/Density-Evolution-AWGN/tree/master/%20Reciprocal-channel%20approximation%20for%20approximation%20of%20Density%20Evolution%20Iterative%20Decoding%20Threshold
contain protograph reciprocal-channel approximation (RCA) for GA like approximation of Density Evolution Iterative Decoding Threshold of  Multi Edge Type (MET) LDPC Codes. Using RCA constructed codes from article  S.-Y. Chung, G.D. Forney, T.J. Richardson, R. Urbanke  On the design of low-density parity-check codes within 0.0045 dB of the Shannon limit, IEEE Communications letters 5 (2), 58-60, 2001. RCA method in detail described at S.-Y.Chung,“On the construction of some capacity-approaching coding schemes,” Ph.D. dissertation, MIT, Cambridge, MA, 2000.  Tested on old (2011) and new version (2018) of Matlab and GNU Octave 5.1 (several times slower than Matlab). RCA-GA is one of the best method for fast DE approximation from runtime and accuracy.
Protograph RCA not require to generate protograph, lift it directly using Simulated Annealing with girth/EMD (https://github.com/Lcrypto/Simulated-annealing-lifting-QC-LDPC) or another QC extension methods. 

2. profgen.m - GA approximation of Density evolution based optimization of degree distribution for LDPC codes under AWGN-channel. Taked from https://www.ece.rice.edu/~debaynas/codes.html Alexandre de Baynast. Due to changes at lsqlin implementation work till R2016(active-set algorithm removed). New version of matlab currently not support. Thank to Dr. Ovinnikov A. for bug report.


3. main.m at subfolder GA - contain Gaussian approximation based Differential Evolution optimization of LDPC Degree distribution from Dr. Chén Zǐ Qiáng (陈紫强), Guilin University of Electronic Technology (桂林电子科技大学，广西桂). 


Very good compare of quality of GA, GA-RCA and Pure DE at paper Sarah J. Johnson, et al "A New Density Evolution Approximation for LDPC and Multi-Edge Type LDPC Codes," in IEEE Transactions on Communications, vol. 64, no. 10, pp. 4044-4056, Oct. 2016, https://arxiv.org/pdf/1605.04665.pdf.

4. Add draft version of Covariance Evolution with Octopus protograph, as generalization of Density Evolution on Finite-Length https://github.com/Lcrypto/Density-Evolution-AWGN/tree/master/Covariance-Evolution. In Dr. Richardson file store tabulated values: 
a,b,VN,CN,Circulant and calculate WER(SNR) = Q(a*(10*log10(SNR)-b)). Using linear interpolation get values for circulant between and calculate WER for them.

5. Quantized Density Evolution from Dr. Andrew W. Eckford   https://github.com/Lcrypto/Density-Evolution-AWGN/tree/master/Quantized%20Density-evolution

6. Very usefull library QuantDMC ver 4 from 2016, Mutual Information optimization for quantized LLRs under Dicrete Memoryless Channel, e.x. non-uniformly quantized AWGN-channel optimization, https://github.com/Lcrypto/Density-Evolution-AWGN/tree/master/QuantDmc-4. Taked from Brian Kurkovskiy  http://www.jaist.ac.jp/is/labs/bits/source. In detail described at B. M. Kurkoski and H. Yagi, “Quantization of binary-input discrete memoryless channels,” IEEE Transactions on Information Theory, vol. 60, no. 8, pp. 4544-4552, August 2014. 
Example of use you can see at https://github.com/Very-Fancy/ldpc-quant
![alt text](https://github.com/Lcrypto/Density-Evolution-AWGN/blob/master/QuantDmc-4/075eng.png)



To get protograph from protograph ensemble (weight profile: weight distribution, rho/lambda) use: 

1. rand_proto.m - generate protograph according to weight distribution of variable and check.
Algorithm probabalistical, can fail. Start script again several time (from 10 to 100K) or use different size of protograph and degree distibution.

2. PEG https://github.com/Lcrypto/classic-PEG-

To lift protograph for certain circulant use:
1. Simulated Annealing with girth/EMD - current state of the art lifting method for optimization of graph properties
https://github.com/Lcrypto/Simulated-annealing-lifting-QC-LDPC
2. To compare you can try Fossorier's approach Guess-and-Test:

https://github.com/Lcrypto/Guess-and-Test-For-CPM-weight-more-1 


and


https://github.com/Lcrypto/-Greedy-Guess-and-Test-method-for-construction-QC-LDPC-codes-with-CPM-of-weigth-more-than-1.

P.S. if you need degree distribution for relay flat fading use [DE for Relay flat-fading](https://github.com/Lcrypto/Density-Evolution-for-relay-flat-fading-channel-)


