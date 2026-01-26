function force_selector = create_force_selector(row, col, mode)
    % row, col: dimensions of the mass grid
    % mode: either 'I' or a vector of indices (e.g., [1, 3])
    
    N = row * col;  % Total number of masses
    total_DOF = 2 * N;  % x and y for each mass
    force_selector = zeros(total_DOF, 1);  % Default: no force anywhere
    
    if ischar(mode) && strcmp(mode, 'I')
        force_selector(:) = 1;  % Full force
    elseif isnumeric(mode)
        % Force only on selected masses
        for i = mode
            if i >= 1 && i <= N
                force_selector(i) = 1;       % x_i
                force_selector(N + i) = 1;   % y_i
            else
                error('Mass index %d is out of bounds. Valid range is 1 to %d.', i, N);
            end
        end
    else
        error('Unsupported mode. Use ''I'' or a numeric array of mass indices.');
    end
end