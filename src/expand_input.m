function u_full = expand_input(u_reduced, force_selector)
    S = actuator_mapping(force_selector);
    assert(numel(u_reduced)==size(S,2),'Reduced input size does not match mask');
    u_full = S*u_reduced(:);
end
