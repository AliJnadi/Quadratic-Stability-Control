function [A, B, b] = gen_sys_matrices(row, col, force_selector, rho)
% gen_sys_matrices will generate A, B matrices for the given row * col
% system of mass-spring-damper.
% Inputs:
% row: number of rows of the grid system.
% col: number of columns of the grid system.
% force_selector: vector represents where forces are applied and by which
% amplitude
% pho: system parameter
% Outputs:
% A system state matrix
% B system control matrix
% b afffine term

    m = rho.m;
    k = rho.k;
    mu = rho.mu;
    lx = rho.lx; 
    ly = rho.ly; 

    % state vector [x, y, dx, dy]
    zero_r_c = zeros(2*row*col);
    eye_r_c = eye(2*row*col);

    A = [zero_r_c, eye_r_c; zeros(2*row*col, 4*row*col)];
    B = [zeros(2*row*col); diag(force_selector)/m];
    b = zeros(4*row*col, 1);

    idx = 1;
    % shifting
    y_shift = row*col;
    d_shift = 2*y_shift;
    
    % Padding with zero for easy dealing with walls;
    r = row + 2;
    c = col + 2; 
    
    for i = 2:r-1
        for j = 2:c-1
            [tempx, bx] = generate_row_x(i, j, r, c, m, k, mu, lx, ly);
            [tempy, by] = generate_row_y(i, j, r, c, m, k, mu, lx, ly);
            
            A(d_shift + idx, :) = tempx;
            A(y_shift + d_shift + idx, :) = tempy;
           
            b(d_shift + idx) = bx;
            b(y_shift + d_shift + idx) = by;
            idx = idx+1;
        end
    end    
end

function [A, b] = generate_row_x(i, j, row, col, m, k, mu, lxx, lyy) 
    % Represent state as matrices for easly deal with walls
    x = zeros(row, col);
    y = zeros(row, col);
    dx = zeros(row, col);
    dy = zeros(row, col);
    
    % forces related to stiffness
    x(i, j-1) = k;
    x(i, j+1) = k;
    
    x(i, j) = -2*k;
    
    % forces related to damping
    dx(i, j-1) = mu;
    dx(i, j+1) = mu;
    
    dx(i, j) = -2*mu; 
    
    %calculate b
    b = calculate_b(row, col, lxx, lyy, x, y, dx, dy);
    
    x = x(2:end-1, 2:end-1); x = reshape(x.',1,[]);
    y = y(2:end-1, 2:end-1); y = reshape(y.',1,[]);
    dx = dx(2:end-1, 2:end-1); dx = reshape(dx.',1,[]);
    dy = dy(2:end-1, 2:end-1); dy = reshape(dy.',1,[]);
    A = [x, y, dx, dy];
    
    A = A / m;
    b = b / m;
end

function [A, b] = generate_row_y(i, j, row, col, m, k, mu, lxx, lyy)  
    x = zeros(row, col);
    y = zeros(row, col);
    dx = zeros(row, col);
    dy = zeros(row, col);
    
    % forces related to stiffness
    y(i-1, j) = k;
    y(i+1, j) = k;
    
    y(i, j) = -2*k;
    
    % forces related to damping
    dy(i-1, j) = mu;
    dy(i+1, j) = mu;
    
    dy(i, j) = -2*mu;
    
    %calculate b
    b = calculate_b(row, col, lxx, lyy, x, y, dx, dy);
    
    x = x(2:end-1, 2:end-1); x = reshape(x.',1,[]);
    y = y(2:end-1, 2:end-1); y = reshape(y.',1,[]);
    dx = dx(2:end-1, 2:end-1); dx = reshape(dx.',1,[]);
    dy = dy(2:end-1, 2:end-1); dy = reshape(dy.',1,[]);
    
    A = [x, y, dx, dy];
    A = A / m;
    
    b = b / m;
end

function b = calculate_b(row, col, lxx, lyy, x, y, dx, dy)
    x = reshape(x.',1,[]);
    y = reshape(y.',1,[]);
    dx = reshape(dx.',1,[]);
    dy = reshape(dy.',1,[]);
    A = [x, y, dx, dy];
    
    [x, y] = meshgrid(0:lxx:(col-1)*lxx, 0:lyy:(row-1)*lyy);
    
    % Nullify interior points
    interior_mask = true(row, col);      % start with all true
    interior_mask(1, :) = false;         % top border
    interior_mask(end, :) = false;       % bottom border
    interior_mask(:, 1) = false;         % left border
    interior_mask(:, end) = false;       % right border
    % Apply mask to nullify inner elements
    x(interior_mask) = 0;
    y(interior_mask) = 0;
    
    dx = zeros(row, col);
    dy = zeros(row, col);
    
    x = reshape(x.',1,[]);
    y = reshape(y.',1,[]);
    dx = reshape(dx.',1,[]);
    dy = reshape(dy.',1,[]);
    s = [x, y, dx, dy]';
    
    b = A*s;
end