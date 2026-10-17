function [S, indices] = actuator_mapping(force_selector)
% Reduced inputs: [x forces on active masses; y forces on active masses].
% Full inputs: [x forces on all masses; y forces on all masses].
assert(isvector(force_selector) && all(ismember(force_selector,[0 1])), ...
    'Actuator mask must be a binary vector');
N = numel(force_selector);
active = find(force_selector(:));
indices = [active; N+active];
S = zeros(2*N,numel(indices));
for j=1:numel(indices)
    S(indices(j),j)=1;
end
end
