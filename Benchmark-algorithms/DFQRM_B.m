function Result = DFQRM_B(fparam, param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%
% This function is an implementation of Algorithm 1 in 
%
% Grapiglia, G.: Worst‑case evaluation complexity of a derivative‑free
% quadratic regularization method. Optimization Letters 18 (2024)
%
% The referred Algorithm 1 it is a derivative-free version of the Gradient 
% Method with Armijo line-search, in which the gradients are approximated 
% by forward finite differences.
%
% Inputs
%
%   fun  :function that computes the objective function values
%   x0   : starting point
%   ep   : target accuracy for the norm of the gradient of the objective
%   nfmax: maximum number of function evaluations allowed
%
% Outputs
%
%   x     : vector yielding minimal function value found;
%   fmin  : function value at x;
%   nf    : number of function evaluations;
%   H     : matrix with two columns such that H(i,2) is the best value
%           of the objective f obtained by the algorithm after 
%           H(i,1) function evaluations.
%
% Author information:
%
% Geovani Nunes Grapiglia
% Université catholique de Louvain, Belgium
% Department of Applied Mathematics
% geovani.grapiglia@uclouvain.be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Initialization
Result = struct;
Result.algname = 'DFQRM';

% Unpack
if isfield(param, 'x0')
    x = param.x0;
    [n m]=size(x);
end
nfmax = param.budget;
ep = 10^-5; % this shouldn't be the reason algorithm terminates.
fun = fparam.f;
[f]=fun(x);
nf=1;
H=[nf f];
sigma=1;
sigma_min=10^(-2);
u=[f];
k=0;
fmin = f;


while (nf<nfmax) && (fmin > 10^(-10))
  k=k+1;
  h=2*ep/(5*sigma*sqrt(n)); 
  [g,H1,nf]=grad_forward(x,fun,f,nf,h);
  H=[H;H1];
  ngrad=norm(g);
  while ngrad<(4/5)*ep
      sigma=2*sigma;
      h=2*ep/(5*sigma*sqrt(n)); 
      [g,H1,nf]=grad_forward(x,fun,f,nf,h);
      H=[H;H1];
      ngrad=norm(g);
  end
  x1=x-(1/sigma)*g;
  nstep1=norm(x1-x);
  c1=f-(sigma/8)*(nstep1^(2));
  [f1]=fun(x1);
  nf=nf+1;
  H=[H;nf f1];
  while f1>c1
    sigma=2*sigma;
    h=2*ep/(5*sigma*sqrt(n)); 
    [g,H1,nf]=grad_forward(x,fun,f,nf,h);
    H=[H;H1];
    ngrad=norm(g);
    while ngrad<(4/5)*ep
      sigma=2*sigma;
      h=2*ep/(5*sigma*sqrt(n)); 
      [g,H1,nf]=grad_forward(x,fun,f,nf,h);
      H=[H;H1];
      ngrad=norm(g);
    end
    x1=x-(1/sigma)*g;
    nstep1=norm(x1-x);
    c1=f-(sigma/8)*(nstep1^(2));
    [f1]=fun(x1);
    nf=nf+1;
    H=[H;nf f1];
  end
  x=x1;
  f=f1;
  u=[u;f];
  sigma=max([sigma/2 sigma_min]);
  fmin = min(u);
  disp(['best f value so far is ', num2str(fmin), ' and step size is ', num2str(1/sigma)])
end
%fmin=min(u);
%plot(u);

% Collate results
Result.objval_seq = H(:,2); %toDo: Fix this so it retains only the best so far.
Result.num_queries = H(:,1);
Result.sol = x;
disp(['Final loss is ', num2str(H(end,2)), ' and the number of queries is ', num2str(H(end,1))])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [g,H,nf]=grad_forward(x,fun,f,nf0,h)

% This function computes an approximation for the gradient of f at x 
% using forward finite differences with stepsize 

H=[];  
[n m]=size(x);
I=eye(n);
g=zeros(n,1);
nf=nf0;
for j=1:n
  z=x+h*I(:,j);
  fz=fun(z);
  nf=nf+1;
  H=[H;nf fz];
  g(j)=(fz-f)/h;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



