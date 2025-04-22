function [sparsity_flag, sparsity_parameter] = check_sparsity(vector, epsilon)
% Function to check approximate sparsity. 
% Daniel McKenzie and Geovani Nunes Grapiglia
% April 2025.
%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
n  = length(vector);
sparsity_flag = false; % default. Will be set to true if sparsity test is passed.
sparsity_parameter = n;


% Sort vector
sorted_vector = sort(abs(vector), 'ascend');
sorted_vector_squared = sorted_vector.^2;

% Compute magnitudes of tail vectors.
partial_ell1_norms = cumsum(sorted_vector);
partial_ell2_norms = cumsum(sorted_vector_squared);

% If the tail, i.e. ||vector - [vector]_s|| is small, in both ell_1 and 
% ell_2 norms, it means that the vector is approximately sparse.
for i = 1:n
    if partial_ell1_norms(i) <= epsilon && partial_ell2_norms(i) <= epsilon
        sparsity_flag = true;
        sparsity_parameter = n - i;
    end
end

end
