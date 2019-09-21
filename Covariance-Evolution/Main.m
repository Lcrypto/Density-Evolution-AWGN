%protograph
protograph=[                      
1 0 0 0 0 1 0 1 0 1 1 0 1 0 0 2 
1 1 0 0 0 0 0 1 1 0 0 1 0 1 1 2 
0 1 1 0 0 0 0 1 0 1 0 1 0 1 0 1 
0 0 1 1 0 0 1 0 1 0 1 0 0 1 1 1  
0 0 0 1 1 0 1 0 0 1 1 0 1 0 1 2 
0 0 0 0 1 1 1 0 1 0 0 1 1 0 0 2  
];

savefile=['Octopus.mat'];    % file to store result
eps=0.2;                 % erasure probability
invstep=4;                        % step^-1 is the step size for Diff Eq.                            
Var=false;                   % flag to compute the variance  
[nc,nv]=size(protograph);
puncture_mask=ones(1,nv); 
puncture_mask(16)=0; %choice symbols for puncturing


input.proto=protograph;
input.puncture_mask=puncture_mask;
input.eps=eps;
input.nv=nv;
input.nc=nc;
input.Var=Var;
input.step=invstep;
input.savefile=savefile;
calc(input)