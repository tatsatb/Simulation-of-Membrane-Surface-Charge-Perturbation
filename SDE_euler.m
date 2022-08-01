function xhat = SDE_euler_deb(bigtheta,PROBLEM,OWNTIME,NUMDEPVARS,NUMSIM,SDETYPE,param,input)
%
% revised 1/22/10 by Yuan:
%   1)  change SEED from 0 to sum(100*clock) as suggested, to get different
%   sequences of pseudopod-random numbers at each run
%   2)  add "param" and "input"


% revised 1/22/10 by Debojyoti:
%   1) Negative values of entities are neglected

% Fixed stepsize Euler-Maruyama scheme for numerical solution of Ito SDE systems.
%
% usage: xhat = SDE_euler(bigtheta)
%
% IN:     bigtheta; complete vector of structural model parameters
%         PROBLEM; the user defined name of the current problem/experiment/example etc. (e.g. 'mySDE')
%         OWNTIME; vector containing the equispaced simulation times sorted in ascending order. 
%                   It has starting simulation-time in first and ending simulation-time in last position. 
%                   Thus OWNTIME(i) - OWNTIME(i-1) = h, where h is the fixed stepsize 
%                   for the numerical intregration (i=2,3,...)
%         NUMDEPVARS; the number of dependent variables, i.e. the SDE dimension
%         NUMSIM; the number of desired simulations for the SDE numerical integration 
%         SDETYPE; the SDE definition: must be 'Ito' 
%         SEED; the seed for the generation of pseudo-random normal variates, i.e. the argument for randn('state',SEED);
%               type 'help randn' for details;
% OUTPUT: xhat; the array of the SDE approximated solution at times OWNTIME 
%
% REFERENCE: [1] Kloeden, Platen and Schurz "Numerical solution of SDE through computer experiments", Springer-Verlag 2nd edition 1997

% Copyright (C) 2007, Umberto Picchini  
% umberto.picchini@biomatematica.it
% http://www.biomatematica.it/Pages/Picchini.html


% input check
%error(nargchk(7, 7, nargin));
narginchk(8, 8);   % Yuan 1/22/10
[n1,n2]=rat(NUMSIM); 
if(~isequal(n2,1) || isempty(NUMSIM) || NUMSIM <= 0)
    error('The number of trajectories NUMSIM must be a positive integer');
end
[n1,n2]=rat(NUMDEPVARS);
if(~isequal(n2,1) || isempty(NUMDEPVARS) || NUMDEPVARS <= 0)
    error('The number of variables NUMDEPVARS must be a positive integer');
end

switch upper(SDETYPE)
    case 'ITO'
        ok = 1;
    otherwise
        error('The Euler-Maruyama scheme is defined only for Ito SDE');
end
        

N = length(OWNTIME);
% handle = waitbar(0,'SDE solver progress...');   % Yuan 1/22/10 %
% commented Debojyoti

%[t, xstart] = feval([PROBLEM, '_sdefile'],OWNTIME(1),[],'init',bigtheta,SDETYPE,NUMDEPVARS,NUMSIM);  % initial conditions
[t, xstart] = feval([PROBLEM, '_sdefile'],OWNTIME(1),[],'init',bigtheta,SDETYPE,NUMDEPVARS,NUMSIM,param,input); % Yuan 1/22/10




XVARS = zeros(N,NUMSIM*NUMDEPVARS);  % the predictions matrix
xstart = xstart';
XVARS(1,:) = xstart([1:size(xstart,1)]' * ones(1,NUMSIM), :)'; % ugly but faster than XVARS(1,:) = repmat(xstart,1,NUMSIM);

% % Control the generation of the pseudo-random standard gaussian draws to
% % get repeatable results
% if(isempty(SEED))
%     randn('state',0);
% else
%     randn('state',SEED);
% end 
% randn('state',sum(100*clock));  % Yuan 1/22/10

rng(100)

for j=2:N
%     waitbar((j-1)/(N-1)); % commented Debojyoti
    % t is inherited as the starting time for this interval
    x = XVARS(j-1, :);           % the value(s) of XVARS at the start of the interval
    h = OWNTIME(j)- t;           % the delta time (end - start) -> fixed size of the step .
    
    % Wiener increments generation with 'antithetic variates' variance reduction method
    Winc = zeros(1,NUMDEPVARS*NUMSIM);
    if(mod(NUMDEPVARS*NUMSIM,2)==0)
        Winc(1 : NUMDEPVARS * NUMSIM/2) = sqrt(h)*randn(1,NUMDEPVARS * NUMSIM/2); % the Wiener increment(s) dWj; 
        Winc(NUMDEPVARS * NUMSIM/2 + 1 : end) = -Winc(1 : NUMDEPVARS * NUMSIM/2); % retrieve the other half of the increments using the 'antithetic variates' variance reduction method;
    else
       % adjustment when (number of simulations * number of variables) is odd
        Winc(1 : round(NUMSIM*NUMDEPVARS/2)) = sqrt(h)*randn(1,round(NUMSIM*NUMDEPVARS/2));  % the Wiener increment(s) dWj; 
        Winc(round(NUMSIM*NUMDEPVARS/2)+1:end) = -Winc(1 : round(NUMSIM*NUMDEPVARS/2)-1);  % retrieve the other half of the increments using the 'antithetic variates' variance reduction method;
    end

    %[f,g] = feval([PROBLEM, '_sdefile'], t, x, [], bigtheta,SDETYPE,NUMDEPVARS,NUMSIM);     % the sdefile output
    [f,g] = feval([PROBLEM, '_sdefile'], t, x, [], bigtheta,SDETYPE,NUMDEPVARS,NUMSIM,param,input);  % Yuan 1/22/10

    XVARS(j , :) = x + f * h + g .* Winc ;  % the Euler-Maruyama scheme for Ito SDEs with DIAGONAL noise
    
    
    ind_neglect = XVARS(j , :)<0;   
    XVARS(j,ind_neglect) = XVARS(j-1,ind_neglect); 
    
    
    t = OWNTIME(j);    % now both t and j refer to the end-of-interval   
end

xhat = XVARS;
% close(handle); % commented Debojyoti
