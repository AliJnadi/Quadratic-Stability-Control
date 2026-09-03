function u_full = expand_input(u_reduced, force_selector)
    N = length(force_selector);       % number of masses
    u_full = zeros(2*N,1);           % full input vector

    % indices of active forces
    idx = find(force_selector);

    % x-components
    u_full(idx) = u_reduced(1:length(idx));

    % y-components
    u_full(N + idx) = u_reduced(length(idx)+1:end);
end