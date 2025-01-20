function hl = data_profile(H,N,starting_value,gate)
%     This subroutine produces a data profile as described in:
%
%     Benchmarking Derivative-Free Optimization Algorithms
%     Jorge J. More' and Stefan M. Wild
%     SIAM J. Optimization, Vol. 20 (1), pp.172-191, 2009.
%
%     The latest version of this subroutine is always available at
%     http://www.mcs.anl.gov/~more/dfo/
%     The authors would appreciate feedback and experiences from numerical
%     studies conducted using this subroutine.
%
%     The subroutine returns a handle to lines in a data profile.
%
%       H contains a three dimensional array of function values.
%         H(f,p,s) = function value # f for problem p and solver s.
%       N is an np-by-1 vector of (positive) budget units. If simplex
%         gradients are desired, then N(p) would be n(p)+1, where n(p) is
%         the number of variables for problem p.
%       gate is a positive constant reflecting the convergence tolerance.
%
%     Argonne National Laboratory
%     Jorge More' and Stefan Wild. January 2008.
%     MODIFIED BY DANIEL MCKENZIE 04/08/2024
%     H is now a two dimensional cell array.
%     H{p,s} = Results struct for problem.

[np,ns] = size(H); % Grab the dimensions

% Produce a suitable history array with sorted entries:
% for j = 1:ns
%     for i = 2:nf
%       H(i,:,j) = min(H(i,:,j),H(i-1,:,j));
%     end
% end

prob_min = inf(np,1);
for i =1:np
    for j=1:ns
        alg_prob_min = min(H{i,j}.objval_seq);
        if alg_prob_min < prob_min(i)
            prob_min(i) = alg_prob_min;
        end
    end
end

prob_max = starting_value;     % The starting value for each problem

% For each problem and solver, determine the number of
% N-function bundles (e.g.- gradients) required to reach the cutoff value
T = zeros(np,ns);
for p = 1:np
    cutoff = prob_min(p) + gate*(prob_max(p) - prob_min(p));
    for s = 1:ns %ns = num solvers/algorithms
        alg_objval_seq = H{p,s}.objval_seq;
        cutoff_point = find(alg_objval_seq <= cutoff,1);
        if (isempty(cutoff_point))
            T(p,s) = NaN;
        else
            alg_num_queries = H{p,s}.num_queries;
            T(p,s) = alg_num_queries(cutoff_point)/N(p);
        end
    end
end

% Other colors, lines, and markers are easily possible:
colors  = ['b' 'r' 'm' 'k' 'c' 'g' 'y'];   lines   = {'-' '-.' '--'};
markers = [ 's' 'o' '^' 'v' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];
labels = cell(ns,1);
for s=1:ns
    labels{s} = H{1,s}.algname;
end

% Replace all NaN's with twice the max_ratio and sort.
max_data = max(max(T));
T(isnan(T)) = 2*max_data;
T = sort(T);

% For each solver, plot stair graphs with markers.
hl = zeros(ns,1);
for s = 1:ns
    [xs,ys] = stairs(T(:,s),(1:np)/np);
    sl = mod(s-1,3) + 1; sc = mod(s-1,7) + 1; sm = mod(s-1,12) + 1;
    option1 = [char(lines(sl)) colors(sc) markers(sm)];    
    hl(s) = plot(xs,ys,option1, 'LineWidth', 2);
    hold on;
end

% Axis properties are set so that failures are not shown, but with the
% max_ratio data points shown. This highlights the "flatline" effect.
legend(labels);
axis([0 1.1*max_data 0 1]);
