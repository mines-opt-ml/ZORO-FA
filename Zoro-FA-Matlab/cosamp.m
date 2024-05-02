function [Sest, err] = cosamp(Phi,u,K,tol,maxiterations)

% Cosamp algorithm
%   Input
%       K : sparsity of Sest
%       Phi : measurement matrix
%       u: measured vector
%       tol : tolerance for approximation between successive solutions. 
%   Output
%       Sest: Solution found by the algorithm
%
% Algorithm as described in "CoSaMP: Iterative signal recovery from 
% incomplete and inaccurate samples" by Deanna Needell and Joel Tropp.
% 


% This implementation was written by David Mary, 
% but modified 20110707 by Bob L. Sturm to make it much clearer,
% and corrected multiple times again and again.
% To begin with, see: http://media.aau.dk/null_space_pursuits/2011/07/ ...
% algorithm-power-hour-compressive-sampling-matching-pursuit-cosamp.html
%
% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Short Disclaimer: this script is for educational purpose only.
% Longer Disclaimer see  http://igorcarron.googlepages.com/disclaimer

% Initialization
Sest = zeros(size(Phi,2),1);
residual = u;
t = 1; 
numericalprecision = 1e-14;
T = [];

while (t <= maxiterations) && (norm(residual)/norm(u) > tol)
  y = abs(Phi'*residual);
  [vals,z] = sort(y,'descend');
  Omega = find(y >= vals(2*K) & y > numericalprecision);
  T = union(Omega,T);
  %b = pinv(Phi(:,T))*u;
  [b,flag] = lsqr(Phi(:,T),u);
  [vals,z] = sort(abs(b),'descend');
  try
    Kgoodindices = find(abs(b) >= vals(K) & abs(b) > numericalprecision);
  catch ME
      T
      % warning('all elements of y are below numerical precision')
      rethrow(ME)
      % Kgoodindices = 1:length(T);
  end
  T = T(Kgoodindices);
  Sest = zeros(size(Phi,2),1);
  b = b(Kgoodindices);
  Sest(T) = b;
  %disp(['norm of solution is', num2str(norm(Sest))])
  residual = u - Phi(:,T)*b;
  %norm_residual = norm(residual);
  %disp(['residual is ', num2str(norm_residual)])
  t = t+1;
  err = norm(residual)/norm(u);
  %disp(['Error is', num2str(err)])
end