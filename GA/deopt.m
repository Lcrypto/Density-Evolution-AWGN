function [FVr_bestmem,S_bestval,I_nfeval,FM_pop] = deopt(S_struct,var_people)
global var_index
%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP    = S_struct.I_NP;
F_weight  = S_struct.F_weight;
F_CR    = S_struct.F_CR;
I_D     = S_struct.I_D;
I_itermax  = S_struct.I_itermax;
I_strategy = S_struct.I_strategy;
iter_GA   = S_struct.iter_GA
%-----Initialize population and some arrays-------------------------------
FM_pop = var_people;last_sigma=0;

FM_popold  = zeros(size(FM_pop)); % toggle population
FVr_bestmem = zeros(1,I_D);% best population member ever
FVr_bestmemit = zeros(1,I_D);% best population member in iteration
I_nfeval   = 0;          % number of function evaluations

%------Evaluate the best member after initialization----------------------

I_best_index = 1;         % start with first population member
S_val(1) =price_evalution(FM_pop(I_best_index,:),iter_GA);
S_bestval = S_val(1);        % best objective function value so far
I_nfeval = I_nfeval + 1;
for k=2:I_NP             % check the remaining members%
 S_val(k) =price_evalution(FM_pop(k,:),iter_GA);
 I_nfeval = I_nfeval + 1;
 if (left_win(S_val(k),S_bestval) == 1)
  I_best_index = k;       % save its location
  S_bestval   = S_val(k);
 end 
end
FVr_bestmemit = FM_pop(I_best_index,:); % best member of current iteration
S_bestvalit = S_bestval;       % best value of current iteration

FVr_bestmem = FVr_bestmemit;      % best member ever

%------DE-Minimization---------------------------------------------
%------FM_popold is the population which has to compete. It is--------
%------static through one iteration. FM_pop is the newly--------------
%------emerging population.----------------------------------------

FM_pm1 = zeros(I_NP,I_D); % initialize population matrix 1
FM_pm2 = zeros(I_NP,I_D); % initialize population matrix 2
FM_pm3 = zeros(I_NP,I_D); % initialize population matrix 3
FM_pm4 = zeros(I_NP,I_D); % initialize population matrix 4
FM_pm5 = zeros(I_NP,I_D); % initialize population matrix 5
FM_bm  = zeros(I_NP,I_D); % initialize FVr_bestmember matrix
FM_ui  = zeros(I_NP,I_D); % intermediate population of perturbed vectors
FM_mui = zeros(I_NP,I_D); % mask for intermediate population
FM_mpo = zeros(I_NP,I_D); % mask for old population
FVr_rot = (0:1:I_NP-1);   % rotating index array (size I_NP)
FVr_rotd = (0:1:I_D-1);   % rotating index array (size I_D)
FVr_rt = zeros(I_NP);        % another rotating index array
FVr_rtd = zeros(I_D);        % rotating index array for exponential crossover
FVr_a1 = zeros(I_NP);        % index array
FVr_a2 = zeros(I_NP);        % index array
FVr_a3 = zeros(I_NP);        % index array
FVr_a4 = zeros(I_NP);        % index array
FVr_a5 = zeros(I_NP);        % index array
FVr_ind = zeros(4);

I_iter = 1;
while ((I_iter < I_itermax))
 FM_popold = FM_pop;         % save the old population
 S_struct.FM_pop = FM_pop;
 S_struct.FVr_bestmem = FVr_bestmem;
 
 FVr_ind = randperm(4);       % index pointer array

 FVr_a1 = randperm(I_NP);         % shuffle locations of vectors
 FVr_rt = rem(FVr_rot+FVr_ind(1),I_NP);  % rotate indices by ind(1) positions
 FVr_a2 = FVr_a1(FVr_rt+1);        % rotate vector locations
 FVr_rt = rem(FVr_rot+FVr_ind(2),I_NP);
 FVr_a3 = FVr_a2(FVr_rt+1);        
 FVr_rt = rem(FVr_rot+FVr_ind(3),I_NP);
 FVr_a4 = FVr_a3(FVr_rt+1);       
 FVr_rt = rem(FVr_rot+FVr_ind(4),I_NP);
 FVr_a5 = FVr_a4(FVr_rt+1);        

 FM_pm1 = FM_popold(FVr_a1,:);      % shuffled population 1
 FM_pm2 = FM_popold(FVr_a2,:);      % shuffled population 2
 FM_pm3 = FM_popold(FVr_a3,:);      % shuffled population 3
 FM_pm4 = FM_popold(FVr_a4,:);      % shuffled population 4
 FM_pm5 = FM_popold(FVr_a5,:);      % shuffled population 5

 for k=1:I_NP               % population filled with the best member
  FM_bm(k,:) = FVr_bestmemit;      % of the last iteration
 end

 FM_mui = rand(I_NP,I_D) < F_CR; % all random numbers < F_CR are 1, 0 otherwise 
 FM_mpo = FM_mui < 0.5;     % inverse mask to FM_mui

 if (I_strategy == 1)              % DE/rand/1
  FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2); % differential variation
  FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;  % crossover
  FM_origin = FM_pm3;
 elseif (I_strategy == 2)            % DE/local-to-best/1
  FM_ui = FM_popold + F_weight*(FM_bm-FM_popold) + F_weight*(FM_pm1 - FM_pm2);
  FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;
  FM_origin = FM_popold;
 elseif (I_strategy == 3)            % DE/best/1 with jitter
  FM_ui = FM_bm + (FM_pm1 - FM_pm2).*((1-0.9999)*rand(I_NP,I_D)+F_weight);       
  FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;
  FM_origin = FM_bm;
 elseif (I_strategy == 4)            % DE/rand/1 with per-vector-dither
  f1 = ((1-F_weight)*rand(I_NP,1)+F_weight);
  for k=1:I_D
    FM_pm5(:,k)=f1;
  end
  FM_ui = FM_pm3 + (FM_pm1 - FM_pm2).*FM_pm5;  % differential variation
  FM_origin = FM_pm3;
  FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;  % crossover
 elseif (I_strategy == 5)             % DE/rand/1 with per-vector-dither
  f1 = ((1-F_weight)*rand+F_weight);
  FM_ui = FM_pm3 + (FM_pm1 - FM_pm2)*f1;    % differential variation
  FM_origin = FM_pm3;
  FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;  % crossover
 else                       % either-or-algorithm
  if (rand < 0.5);               % Pmu = 0.5
    FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);% differential variation
    FM_origin = FM_pm3;
  else                     % use F-K-Rule: K = 0.5(F+1)
    FM_ui = FM_pm3 + 0.5*(F_weight+1.0)*(FM_pm1 + FM_pm2 - 2*FM_pm3);
  end
  FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;  % crossover  
 end
 
%-----Optional parent+child %selection-----------------------------------------
%-----Select which vectors are allowed to enter the new population------------
 for k=1:I_NP
   FM_ui(k,:)=abs(FM_ui(k,:));tmp=sum(FM_ui(k,:));           
 if tmp<1
   S_tempval= price_evalution2(FM_ui(k,:),iter_GA,S_val(k));
   I_nfeval = I_nfeval + 1;   
  if S_tempval>S_val(k) 
    FM_pop(k,:) = FM_ui(k,:);          % replace old vector with new one (for new iteration)
    S_val(k) = S_tempval;           % save value in "cost array"
    %----we update S_bestval only in case of success to save time-----------
    if S_tempval>S_bestval 
      S_bestval = S_tempval;          % new best value
      FVr_bestmem = FM_ui(k,:);        % new best parameter vector ever
    end
   end
 end %if
 end % for k = 1:NP

 FVr_bestmemit = FVr_bestmem;   % freeze the best member of this iteration for the coming
 if S_bestval>last_sigma
 I_iter = I_iter + 1
 s=sprintf('%dx%d_iter_sigma.txt',I_NP,I_D);
 fid = fopen(s,'a');
 fprintf(fid,'\n');
 fprintf(fid,'%s %6.3E\n','iter =',I_iter);
 fprintf(fid,'%s %6.8E\n','sigma =',S_bestval); 
 fclose(fid);
 save rate09 S_bestval FVr_bestmem FM_pop var_index
 end
last_sigma=S_bestval;
end %---end while ((I_iter < I_itermax) ...
