function [H,exactRate,nk,ERROR_FLAG] = hgen(R,N,Lambda,Deg_Lambda,Rho,Deg_Rho,seed,DEBUG_FLAG,FileName_sparse);
% [H_SPARSE,EXACT_RATE,NK,ERROR_FLAG] = HGEN(R,N,LAMBDA,DEG_LAMBDA,RHO,DEG_RHO,SEED,DEBUG_FLAG, FILENAME_SPARSE)
% generates a random Parity-Check Matrix of rate R and length N for a given
% profile
% LAMBDA = [LAMBDA(1),...,LAMBDA(.)] is the polynomial representing the proportion of edges LAMBDA(i) connected to
% variables nodes of degree DEG_LAMBDA(i)
% DEG_LAMBDA = [2,3,...] is the vector representing all degrees of the variables-nodes in the code
% RHO = [RHO(1),...,RHO(.)] is the polynomial representing the proportion of edges RHO(i) connected to
% check nodes of degree DEG_RHO(i)
% DEG_RHO = [2,3,...] is the vector representing all degrees of the check-nodes in the code
% SEED seed number to initialize the random generator
% DEBUG_FLAG 0 to desactive it, 1 to print some crucial variables
% FILENAME to use the default name (gallager.n.nk.rate.dat), use []
% EXACT_RATE the exact rate given by the parity-check matrix (as close as
% possible of R)
% NK dimension of the parity-check equation H (nk x n)
% EXAMPLE:
%R = .6;N = 6500;Lambda = [0.2096540000 0.4024980000 0.3878480000];Deg_Lambda = [2 3 10];Rho = [1];Deg_Rho = [9];seed = 0;DEBUG_FLAG = 1;FileName = [];

ERROR_FLAG = 0;
FALSE_ = 0;
TRUE_ = 1;
rand('state',seed);

if (R-floor(100*R)/100<2e-2);
    R = floor(100*R)/100;
end;
CodeNumber = 0;
Int_0_1_Lambda = sum(Lambda./Deg_Lambda);
Int_0_1_Rho = sum(Rho./Deg_Rho);
M = round(N*(1-round(1e16*(1-Int_0_1_Rho/Int_0_1_Lambda))/1e16));
Proportion_Lambda = (Int_0_1_Lambda)*Lambda./Deg_Lambda/sum((Int_0_1_Lambda)*Lambda./Deg_Lambda);
Proportion_Rho = (Int_0_1_Rho)*Rho./Deg_Rho/sum((Int_0_1_Rho)*Rho./Deg_Rho);
Nb_Variable_Nodes = round((N/Int_0_1_Lambda)*Lambda./Deg_Lambda);
Nb_Check_Nodes = round((M/Int_0_1_Rho)*Rho./Deg_Rho);

% First step: add nb_row_toaddtoright rows of the lowest connection degree;
nb_row_toaddtoright = M-sum(Nb_Check_Nodes);
Nb_Check_Nodes(1) = Nb_Check_Nodes(1)+nb_row_toaddtoright;
% Update the total number of ones on the right side
nb1_right = sum(Nb_Check_Nodes.*Deg_Rho);
nb1_left = sum(Nb_Variable_Nodes.*Deg_Lambda);
nb_col_toaddtoleft = N-sum(Nb_Variable_Nodes);
if length(Lambda)==1 & nb_col_toaddtoleft~=0;
    disp('Cannot generate this matrix');
    return;
else;
    if (nb1_right ~= nb1_left | nb_col_toaddtoleft ~= 0);
        x = [Deg_Lambda(1),Deg_Lambda(2);1,1]\[nb1_right-nb1_left;nb_col_toaddtoleft];
        Nb_Variable_Nodes(1) = Nb_Variable_Nodes(1)+x(1);
        Nb_Variable_Nodes(2) = Nb_Variable_Nodes(2)+x(2);    
    end;
end;

exactRate = (sum(Nb_Variable_Nodes)-sum(Nb_Check_Nodes))/sum(Nb_Variable_Nodes)
KK = sum(Nb_Variable_Nodes)-sum(Nb_Check_Nodes);
nk = M;

if sum(Nb_Variable_Nodes.*Deg_Lambda) ~= sum(Nb_Check_Nodes.*Deg_Rho);
    disp('Cannot generate this matrix');
    return;
end;



if DEBUG_FLAG;
    String = [];
	for l = 1:length(Deg_Lambda);
		String = [String,'+',num2str(Deg_Lambda(l)),'*',num2str(Nb_Variable_Nodes(l))];
	end;
	String = [String,' (',num2str(sum(Nb_Variable_Nodes.*Deg_Lambda)),')','~='];
	for l = 1:length(Deg_Rho);
		String = [String,'+',num2str(Deg_Rho(l)),'*',num2str(Nb_Check_Nodes(l))];
	end;
	String = [String,' (',num2str(sum(Nb_Check_Nodes.*Deg_Rho)),')'];
	disp(String);
	disp(['Add ',num2str(N-sum(Nb_Variable_Nodes)),' columns and ',num2str(M-sum(Nb_Check_Nodes)),' rows']);
	
	disp('Press any key to continue and adjust the weights by hands if necessary');pause;
    %end;
end;

CumSum_Nb_Variable_Nodes = cumsum(Nb_Variable_Nodes);
Nb_Elements_Per_Column = zeros(1,N);
Nb_Elements_Per_Column(1:CumSum_Nb_Variable_Nodes(1)) = Deg_Lambda(1);
for l = 1:length(Deg_Lambda)-1;
    Nb_Elements_Per_Column(CumSum_Nb_Variable_Nodes(l)+1:CumSum_Nb_Variable_Nodes(l+1)) = Deg_Lambda(l+1);
end;

CumSum_Nb_Check_Nodes = cumsum(Nb_Check_Nodes);
Nb_Elements_Per_Row = zeros(1,M);
Nb_Elements_Per_Row(1:CumSum_Nb_Check_Nodes(1)) = Deg_Rho(1);
for l = 1:length(Deg_Rho)-1;
	Nb_Elements_Per_Row(CumSum_Nb_Check_Nodes(l)+1:CumSum_Nb_Check_Nodes(l+1)) = Deg_Rho(l+1);
end;
Random_Nb_Elements_Per_Row = Nb_Elements_Per_Row(randperm(length(Nb_Elements_Per_Row)));


% disp('---------------------- Parity-check matrix (Hr) generation... ------------------');

CYCLESREMOVAL_ = FALSE_;

while ~CYCLESREMOVAL_;

Max_Deg_Rho = max(Deg_Rho);
Index_Temp = 0;
Vec_Temp = zeros(sum(Nb_Check_Nodes.*Deg_Rho),1);
Vec_Temp(Index_Temp+1:Index_Temp+Deg_Lambda(1)*Nb_Variable_Nodes(1)) = repmat([1:CumSum_Nb_Variable_Nodes(1)],1,Deg_Lambda(1));
Index_Temp = Index_Temp+Deg_Lambda(1)*Nb_Variable_Nodes(1);
for l=2:length(Deg_Lambda);
	Vec_Temp(Index_Temp+1:Index_Temp+Deg_Lambda(l)*Nb_Variable_Nodes(l)) = repmat([CumSum_Nb_Variable_Nodes(l-1)+1:CumSum_Nb_Variable_Nodes(l)],1,Deg_Lambda(l));
	Index_Temp = Index_Temp+Deg_Lambda(l)*Nb_Variable_Nodes(l);
end;
Random_Columns_Index = Vec_Temp(randperm(length(Vec_Temp)));

Hr = zeros(M,max(Deg_Rho));
Index = 1;
for n = 1:M;
	Hr(n,1:Random_Nb_Elements_Per_Row(n)) = Random_Columns_Index(Index:Index+Random_Nb_Elements_Per_Row(n)-1).';
	Index = Index+Random_Nb_Elements_Per_Row(n);
end;
%Hr = reshape(Vec_Temp(randperm(length(Vec_Temp))),M,Deg_Rho);

% disp('---------- Parity-check matrix (Hr) generation done ------------');

% ---------- Removal of cycles of length 2... --------------------

% Validation
for n = 1:M;
    if DEBUG_FLAG;
        if mod(n,floor(M/100))==0;disp([num2str(round(n/M*100)),'%']);end;
    end;
    Hrn = Hr(n,1:Random_Nb_Elements_Per_Row(n));
    [B,In,Jn] = unique(Hrn);
	IndicesToSwitch = setdiff([1:Random_Nb_Elements_Per_Row(n)],In);
    NbIndicesToSwitch = length(IndicesToSwitch);
    if n+NbIndicesToSwitch<=M;
        for d = 1:NbIndicesToSwitch;
            IndexToSwitch = IndicesToSwitch(d);
            IndexToSwitch2 = floor(Random_Nb_Elements_Per_Row(n+d)*rand)+1;
            tmp = Hr(n,IndexToSwitch);
            Hr(n,IndexToSwitch) = Hr(n+d,IndexToSwitch2);
            Hr(n+d,IndexToSwitch2) = tmp;
        end;
    else;
        for d = 1:NbIndicesToSwitch;
            IndexToSwitch = IndicesToSwitch(d);
            IndexToSwitch2 = floor(Random_Nb_Elements_Per_Row(n-d)*rand)+1;
            tmp = Hr(n,IndexToSwitch);
            Hr(n,IndexToSwitch) = Hr(n-d,IndexToSwitch2);
            Hr(n-d,IndexToSwitch2) = tmp;
        end;
    end;
    
% 	if Diff ~= 0;
% 		switch Diff;
%             case 1;
% 			    %disp('Attempt to remove cycle(s) of length 2');
% 			    [B,In,Jn] = unique(Hr(n,1:Random_Nb_Elements_Per_Row(n)));
%                 
% 			    for l = 1:Random_Nb_Elements_Per_Row(n);
% 				    if length(find(In==l))==0;
% 					    l_eliminate = l;
% 				    end;
% 			    end;
%                 if ~l_eliminate;
%                     disp('FAILURE --- Matrix is too full -- Please restart the program');
%                     CYCLESREMOVAL_ = FALSE_
%                 end;
% 			    temp = Hr(n,l_eliminate);
%                 if n<M;
% 			        l_switch = floor(Random_Nb_Elements_Per_Row(n+1)*rand)+1;
%     			    Hr(n,l_eliminate) = Hr(n+1,l_switch);
% 			        Hr(n+1,l_switch) = temp;
%                 else;
%                     l_switch = floor(Random_Nb_Elements_Per_Row(n-1)*rand)+1;
%     			    Hr(n,l_eliminate) = Hr(n-1,l_switch);
% 			        Hr(n-1,l_switch) = temp;
%                 end;
%             case {2,3,4,5};
%                 [B,In,Jn] = unique(Hr(n,1:Random_Nb_Elements_Per_Row(n)));
%                 l_eliminate = [];
% 			    for l = 1:Deg_Rho;
% 				    if length(find(In==l))==0;
% 					    l_eliminate = [l_eliminate,l];
% 				    end;
% 			    end;
%                 if length(l_eliminate) < Diff;
%                     disp('FAILURE --- Matrix is too full -- Please restart the program');
%                     CYCLESREMOVAL_ = FALSE_
%                 end;
%                 
%             otherwise;
%                 CYCLESREMOVAL_ = FALSE_
%                 Diff
% 			    %disp('Bad draw. Restart the program');
% 		end;
end;

%disp('Second pass for checking removal of cycle of length 2');
counter = 0;
for n = 1:M;
    %if DEBUG_FLAG;
        if mod(n,floor(M/100))==0;disp([num2str(round(n/M*100)),'%']);end;
    %end;
	Diff = length(Hr(n,1:Random_Nb_Elements_Per_Row(n)))-length(unique(Hr(n,1:Random_Nb_Elements_Per_Row(n))));
	if Diff ~= 0;
		counter = counter+1;
        CYCLESREMOVAL_ = FALSE_;
		disp('Cycle of length 2 still here. The program will automatically restart');pause(2);
	end;
end;
if counter==0;
	disp('Succeed. No more cycles of length 2.')
    CYCLESREMOVAL_ = TRUE_;
end;

end; % end of while ~CYCLESREMOVAL_;

%disp('---------------------- Removal of cycles of length 2 done-------------------------');

if DEBUG_FLAG;
disp('---------------------- Check out the profile of Hr... ----------------------------');

profilec = zeros(1,N);
for n = 1:N;
	if mod(n*10,N) == 0;disp([num2str(n*100/N),'%']);end; 
	[I,J] = find(Hr==n);
	profilec(n) = length(I);
end;
figure(1);plot(profilec,'k-');hold on;

disp('---------------------- Check out the profile of Hr done---------------------------');

end; % end of DEBUG_FLAG;

%disp('---------------------- Generation of Hc... ---------------------------------------');


Hc = zeros(N,max(Deg_Lambda));
for n = 1:N;
    %if DEBUG_FLAG; 
        if mod(n,floor(N/100))==0;disp([num2str(round(n/N*100)),'%']);end;
    %end;
	[I,J] = find(Hr==n);
	Hc(n,1:length(I))=I.';
end;

%disp('---------------------- Generation of Hc done --------------------------------------');



if DEBUG_FLAG;
disp('---------------------- Check out the profile of Hc... ----------------------------');
profiler = zeros(1,N);
for n = 1:M;
	if mod(n*10,M) == 0;disp([num2str(n*100/(M)),'%']);end; 
	[I,J] = find(Hc==n);
	profiler(n) = length(I);
end;
plot(profiler,'r-');hold on;
disp('---------------------- Check out the profile of Hc done---------------------------');
end;% end of DEBUG_FLAG;

if DEBUG_FLAG;
disp('---------------------- Check the rank of H...--------------------------------------');
if N<=1152;
    Hnonsparse = zeros(M,N);
    for ink = 1:M;
        Hnonsparse(ink,Hr(ink,:)) = 1;
    end;
    Hnonsparsegf = gf(Hnonsparse);
    disp('H gf done');
    rankH = rank(Hnonsparsegf);
    disp(['rank(H,',num2str(length(Hnonsparse(:,1))),'x',num2str(length(Hnonsparse(1,:))),') = ',num2str(rankH)]);
else;
    disp('Dimensions of the PCM are too large');
end;
disp('---------------------- Check the rank of H done--------------------------------------');
end;


%disp('---------------------- Printing H in file (matlab sparse format)...--------------------------------------');
if length(FileName_sparse) == 0;
    FileName_sparse = ['gallager.sparse.',num2str(N),'.',num2str(KK),'.',num2str(CodeNumber),'.mat'];
end;
ii = zeros(1,sum(Nb_Elements_Per_Column));
jj = ii;
kkk = 1;
for iii=1:N;
    %if DEBUG_FLAG; 
            if mod(iii,floor(N/100))==0;disp([num2str(round(iii/N*100)),'%']);end;
    %end;
      for jjj=1:Nb_Elements_Per_Column(iii);
      jj(kkk) = iii;
      ii(kkk) = Hc(iii,jjj);
      ss = 1;
      kkk = kkk+1 ; 
   end
end
H = sparse(ii,jj,ss,M,N);
save(FileName_sparse,'H');

%disp('---------------------- Printing H in file done -------------------------------------');

