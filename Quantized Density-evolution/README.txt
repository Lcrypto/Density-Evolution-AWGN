Density Evolution version 0.1.1
Copyright (C) 2003 by Andrew W. Eckford

***Licensing information is contained in the file license.txt
***You must read and agree to the terms of the license before 
***using this software

Revision History

Version 0.1 - January 17, 2003 - First release
Version 0.1.1 - March 14, 2003 - Bug fix in new_chk_overflow.m

This package contains MATLAB scripts which implement Richardson and Urbanke's 
density evolution technique to find the ultimate performance of LDPC codes
in memoryless channels.

Use:

result_pe = de_regular(chan,iter,ext,mapping,stop_pe,dv,dc);

where:

result_pe : is the result vector of length iter, giving the
probability of error after each iteration.

chan : is a vector containing the PDF of the channel message, which
corresponds to the requirements of ext (see below).  Note that 
since chan is a PDF, sum(chan)*increment = 1.

iter : is the maximum number of iterations

ext : is a vector of the form [base increment length], describing the
quantization of the channel and extrinsic PDFs, where "base" is the
least quantization step, "increment" is the difference between adjacent
quantization steps, and "length" is the number of quantization steps.
For example, [-30 0.01 6001] represents a quantization of
-30, -29.99, -29.98, ..., 29.99, 30.

mapping : is a vector of the form [base increment length], used 
internally to represent the density of the log-magnitude of the tanh
of a check message.  These values are always negative since |tanh(x)|<1.
A value I have often used is [-10 0.0002 50000].

stop_pe : is a probability.  If the probability of error goes below stop_pe,
the routine exits, rather than performing the remaining operations.  Set to 
zero to disable.

dv : is the variable degree.

dc : is the check degree.

Example (five iterations for a BSC with inversion probability 0.08394):

>> chan = zeros(1,6001);
>> chan(3240) = 91.606;
>> chan(2762) = 8.394;
>> ext = [-30 0.01 6001];
>> mapping = [-10 0.0002 50000];
>> dv = 3;
>> dc = 6;
>> iter = 5;
>> stop_pe = 1e-5;
>> result_pe = de_regular(chan,iter,ext,mapping,stop_pe,dv,dc)

result_pe =

    0.0839    0.0788    0.0756    0.0730    0.0710



For irregular LDPC just use irregular.m

vard(1,:)=[0 0.2895    0.3158 0 0  0.3947];
chkd(1,:)=[ 0 0 0 0 0 0.9032  0.0968];

