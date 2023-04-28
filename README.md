# Density Evolution, Covariance Evolution for AWGN Symmetrical and Assymetrical channels Matlab and Python implementation
1. The "main.m" file in the subfolder https://github.com/Lcrypto/Density-Evolution-AWGN/tree/master/ Reciprocal-channel approximation for approximation of Density Evolution Iterative Decoding Threshold contains a protograph reciprocal-channel approximation (RCA) for GA-like approximation of Density Evolution Iterative Decoding Threshold of Multi Edge Type (MET) LDPC Codes. The RCA method constructs codes using the method outlined in S.-Y. Chung, G.D. Forney, T.J. Richardson, R. Urbanke's article titled "On the design of low-density parity-check codes within 0.0045 dB of the Shannon limit," which was published in IEEE Communications Letters in 2001. The RCA method is thoroughly described in S.-Y.Chung's Ph.D. dissertation titled "On the construction of some capacity-approaching coding schemes" from MIT in 2000.

The code has been tested on both old (2011) and new (2018) versions of Matlab as well as GNU Octave 5.1 (which is several times slower than Matlab). RCA-GA is one of the best methods for fast DE approximation from runtime and accuracy.

Protograph RCA does not require the generation of protographs; rather, it can be lifted directly using Simulated Annealing with girth/EMD (https://github.com/Lcrypto/Simulated-annealing-lifting-QC-LDPC) or another QC extension method. The "main.m" file provides an example of threshold estimation using the Octopus protograph.

![alt text](https://github.com/Lcrypto/Density-Evolution-AWGN/blob/master/Octopus.png)
 
 
2. "profgen.m" provides GA-approximation of Density Evolution based optimization of degree distribution for LDPC codes under AWGN-channel. The code was taken from Alexandre de Baynast's website at https://www.ece.rice.edu/~debaynas/codes.html. Due to changes in the lsqlin implementation, the code works only until R2016 (active-set algorithm removed), and the new version of Matlab is currently not supported. We would like to thank Dr. Ovinnikov A. for reporting the bug.


3. The "main.m" file in the GA subfolder contains Gaussian approximation-based Differential Evolution optimization of LDPC degree distribution by Dr. Chén Zǐ Qiáng (陈紫强) from Guilin University of Electronic Technology (桂林电子科技大学，广西桂).


The paper titled "A New Density Evolution Approximation for LDPC and Multi-Edge Type LDPC Codes" by Sarah J. Johnson et al., published in IEEE Transactions on Communications in October 2016 (https://arxiv.org/pdf/1605.04665.pdf), provides an excellent comparison of the quality of GA, GA-RCA, and pure DE methods.

4. A draft version of Covariance Evolution, which provides a solution to the peeling decoder differential equation with Octopus protograph and is a generalization of Density Evolution on Finite-Length, is now available at https://github.com/Lcrypto/Density-Evolution-AWGN/tree/master/Covariance-Evolution. In Dr. Richardson's file, tabulated values such as a, b, VN, CN, and Circulant are stored. WER(SNR) = Q(a*(10*log10(SNR)-b)) is calculated using these values. Linear interpolation is used to obtain values for circulant between and to calculate WER for them.

5. Dr. Fan Zhang's paper titled "Covariance Evolution Estimation of Irregular LDPC List Sum-Product Decoder under Q-ary Symmetric Channel" (https://arxiv.org/abs/0806.3243) provides an estimate of irregular LDPC List Sum-Product decoder under q-ary symmetric channel. The MATLAB source code for this estimation is available at https://github.com/Lcrypto/Density-Evolution-AWGN/tree/master/Covariance_Evolution_of_LIST_BP_under_q-ary_channel.

6. Dr. Andrew W. Eckford has developed Quantized Density Evolution, which is available at https://github.com/Lcrypto/Density-Evolution-AWGN/tree/master/Quantized Density-evolution.

7. The QuantDMC ver 4 library, developed in 2016, provides Mutual Information optimization for quantized LLRs under Discrete Memoryless Channel. For example, it includes non-uniformly quantized AWGN-channel optimization, which can be found at https://github.com/Lcrypto/Density-Evolution-AWGN/tree/master/QuantDmc-4. Brian Kurkovskiy's website at http://www.jaist.ac.jp/is/labs/bits/source hosts the code. The paper titled "Quantization of binary-input discrete memoryless channels" by B. M. Kurkoski and H. Yagi, published in IEEE Transactions on Information Theory in August 2014, provides an in-depth description of the library.

An example of how to use the QuantDMC ver 4 library can be found at https://github.com/Very-Fancy/ldpc-quant. An image from the library is shown below: 
![alt text](https://github.com/Lcrypto/Density-Evolution-AWGN/blob/master/QuantDmc-4/075eng.png)



8. The Density Evolution for Asymmetric Memoryless Channels has been implemented in Python and can be found at https://github.com/Lcrypto/Density-Evolution-AWGN/tree/master/density_evolution_ga_python. Differential evolution is used to find optimal solutions constrained by maximum allowed degree distributions. The Gaussian Approximation implementation is used to construct Multi-edge or Protograph QC-LDPC codes for Richardson-Urbanke Flarion-Qualcomm and Divsalar Jet Propulsion Lab, respectively. It is important to note that the terms "Multi-edge" and "Protograph" refer to the same QC-LDPC codes.

9. The "Reported Thresholds and BER Performance for LDPC and LDPC-Like Codes" by Sarah J. Johnson provides valuable information on Density Evolution threshold. The report can be found at https://github.com/Lcrypto/Density-Evolution-AWGN/blob/master/Reported Thresholds and BER for LDPC.pdf.




To obtain a protograph from a protograph ensemble (weight profile: weight distribution, rho/lambda), you can use various methods.

"rand_proto.m" generates a protograph according to the weight distribution of the variable and check. The algorithm is probabilistic and may fail. To improve success rate, you can try running the script multiple times (from 10 to 100K) or using different sizes of protographs and degree distributions.

PEG C++ (https://github.com/Lcrypto/classic-PEG-/tree/master/classic_PEG), PEG+ACE Matlab (https://github.com/Lcrypto/classic-PEG-/blob/master/ProgressiveEdgeGrowthACE.m), and QC PEG with ACE C++ (https://github.com/Lcrypto/classic-PEG-/tree/master/QC-LDPC ACE-PEG) are other options for obtaining protographs.

To lift a protograph for a certain circulant, you can use:

Simulated Annealing with girth/EMD - this is the current state-of-the-art lifting method for optimizing graph properties (https://github.com/Lcrypto/Simulated-annealing-lifting-QC-LDPC).

Quasi-Cyclic PEG with ACE (https://github.com/Lcrypto/classic-PEG-/tree/master/QC-LDPC ACE-PEG) is a previous best solution.

Fossorier's approach Guess-and-Test can also be used for comparison:

https://github.com/Lcrypto/Guess-and-Test-For-CPM-weight-more-1 


and


https://github.com/Lcrypto/-Greedy-Guess-and-Test-method-for-construction-QC-LDPC-codes-with-CPM-of-weigth-more-than-1.

P.S. Note that if you need the degree distribution for relay flat fading, you can use DE for Relay flat-fading which provides the necessary information. [DE for Relay flat-fading](https://github.com/Lcrypto/Density-Evolution-for-relay-flat-fading-channel-)


