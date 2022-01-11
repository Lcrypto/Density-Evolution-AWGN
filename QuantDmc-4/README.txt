QuantDMC - Matlab library for quantizing discrete memoryless channels.

This software contains implementations of the algorithms for quantizing 
discrete memoryless channels, where mutual information is the objective 
function. The main library features are Matlab implementations of:

* Dynamic programming algorithm, for finding optimal quantizers for 
  binary-input DMCs [1],
* KL-means algorithm, for efficiently finding good quantizers for arbitrary 
  DMCs [2],
* Greedy combining algorithm, an earlier algorithm for arbitrary DMCs, 
  which can sometimes outperform KL-means algorithm at the expense of 
  higher complexity [3]. 

If you use this software and write a paper, please cite relevant publications:

[1] B. M. Kurkoski and H. Yagi, "Quantization of binary-input discrete 
    memoryless channels," IEEE Transactions on Information Theory, vol. 60, 
    no. 8, pp. 4544-4552, August 2014. http://dx.doi.org/10.1109/TIT.2014.2327016

[2] A. Zhang and B. M. Kurkoski, "KL Means Algorithm for Quantization of 
    Discrete Memroyless Channels," Proceedings of International Symposium on 
    Information Theory, July 2016.

[3] B. M. Kurkoski, K. Yamaguchi and K. Kobayashi, "Noise Thresholds for 
    Discrete LDPC Decoding Mappings," Proceedings of IEEE Global 
    Communications Conference, December 2008.

This software is provided under an MIT license to supplement papers on 
quantization of discrete memoryless channels. See LICENSE for details.


INSTALLATION 

Download the zip file from
   http://www.kurkoski.org/source

To install, extract the tar file:
   % unzip QuantDmc-N.zip
where N is the version number.

In Matlab, cd to the directory and STARTUP modifies your search path:
   >> cd QuantDmc-N
   >> startup
   >> help quantDmc

The library includes a MEX function and binaries for some platforms.
If your platform is not included, compile using:
   >> mexify

To run an example
   >> cd example
   >> biAwgnQuantization

Get help:
   >> help quantDmc


VERSION HISTORY

Version 4, 2016 March 11

    * Added QUANTDMCKLMEANS function to implement the KL means clustering 
      algorithm for quantizing non-binary input DMCs.
    * Added QUANTDMCGREEDY function to implement the greedy combining 
      algorithm for quantizing non-binary input DMCs
    * Modified CHANNELSORT to accept non-normalized input, so that 
      ChannelSort(rand(2,M)) produces a random channel.
    * Changed filenames.  Functions have leading lowercase.  Now uses 
      "Bi" in filename to identify binary-input quantization algorithms.
    * Added examples using QUANTDMCKLMEANS, QUANTDMCGREEDY
    * Removed LDPC decoding functions.


Version 3, 2012 June 21 
    
    * Added functions QdeConvergence and QdeThreshold to find noise 
      threshold for quantized messages LDPC decoder using maximization of
      mutual information.


Version 2, 2012 May 31
    
    * Added BiAwgn2Dmc and ExBiAwgn
    * The mexify function is more robust to path changes
    * Fixed a bug that prevented compiling with Microsoft Visual C++ (use 
      DBL_MAX/DBL_MIN from float.h instead of a fixed constant and use 
      log(x) / log(2) instead of log2(x) ).


Version 1, 2012 May 23

    * Initial release.  Contains Matlab and C/MEX functions to find the 
      optimal quantizer for a DMC in the sense of maximizing mutual 
      information.

