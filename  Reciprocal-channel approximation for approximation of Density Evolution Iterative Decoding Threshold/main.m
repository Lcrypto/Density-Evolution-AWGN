Protograph=[
2	1	0	0	0	0	1	0	1	0	1	1	0	1	0	0
2	1	1	0	0	0	0	0	1	1	0	0	1	0	1	1
1	0	1	1	0	0	0	0	1	0	1	0	1	0	1	0
1	0	0	1	1	0	0	1	0	1	0	1	0	0	1	1
2	0	0	0	1	1	0	1	0	0	1	1	0	1	0	1
2	0	0	0	0	1	1	1	0	1	0	0	1	1	0	0
]; % example of Octopus MET QC-LDPC Codes R=2/3, 0.1 dB gain compare AR4JA codes
[C V]=size(Protograph);
iterations=250; %iterations numbers
Snr_start = 2; % starting search point snr value
Punctured_VN =1; %Number of Punctured Nodes from 1st columns in protograph
InfVNs=V-C; %Number of information variable nodes which consider in approximation
%if it some cascade construction (check nodes important) should be set by hand
Rate=(V-C)/(V-Punctured_VN) %Rate of code only for EB_No_result
M=2; % Modulation 2 bpsk, 4 qpsk
snr_result = RCA_threshold(Protograph,Snr_start,Punctured_VN,InfVNs,iterations)
EB_No_result= snr_result-10*log10(log2(M)*Rate)
