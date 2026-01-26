function f_x = non_linear_sys(state, row, col, rho)
    m = rho.m;

    kxy = rho.kxy;
    mu_xy = rho.muxy;
    lx = rho.lx; 
    ly = rho.ly; 
    lxy = sqrt(rho.lx^2 + rho.ly^2);

    dstate_x = zeros(1, row+col);
    dstate_y = zeros(1, row+col);
    ddstate_x = zeros(1, row+col);
    ddstate_y = zeros(1, row+col);
    
    % Represents the main nodes of the system     
    [full_x, full_y] = meshgrid(0:lx:(col + 1)*lx, 0:ly:(row + 1)*ly);
    full_dx = zeros(row+2, col+2);
    full_dy = zeros(row+2, col+2);
    
    % Extract system state     
    [x, y, dx, dy] = state_extractor(state);
    
    x = reshape(x, row, col)';
    y = reshape(y, row, col)';
    dx = reshape(dx, row, col)';
    dy = reshape(dy, row, col)';
    
    full_x(2:end-1, 2:end-1) = x;
    full_y(2:end-1, 2:end-1) = y;
    full_dx(2:end-1, 2:end-1) = dx;
    full_dy(2:end-1, 2:end-1) = dy;
    
    idx = 1;
    for i = 2 : row + 1
        for j = 2 : col + 1
            r1 = sqrt((full_x(i, j) - full_x(i-1, j-1))^2 + (full_y(i, j) - full_y(i-1, j-1))^2);
            r2 = sqrt((full_x(i, j) - full_x(i-1, j+1))^2 + (full_y(i, j) - full_y(i-1, j+1))^2);
            r3 = sqrt((full_x(i, j) - full_x(i+1, j-1))^2 + (full_y(i, j) - full_y(i+1, j-1))^2);
            r4 = sqrt((full_x(i, j) - full_x(i+1, j+1))^2 + (full_y(i, j) - full_y(i+1, j+1))^2); 
            
%             dr1 = sqrt((full_dx(i, j) - full_dx(i-1, j-1))^2 + (full_dy(i, j) - full_dy(i-1, j-1))^2); 
%             dr2 = sqrt((full_dx(i, j) - full_dx(i-1, j+1))^2 + (full_dy(i, j) - full_dy(i-1, j+1))^2);
%             dr3 = sqrt((full_dx(i, j) - full_dx(i+1, j-1))^2 + (full_dy(i, j) - full_dy(i+1, j-1))^2); 
%             dr4 = sqrt((full_dx(i, j) - full_dx(i+1, j+1))^2 + (full_dy(i, j) - full_dy(i+1, j+1))^2);

            f_x = 0;
            % On x direction 
            f_x = f_x - kxy*(r1 - lxy)* (full_x(i, j) - full_x(i - 1, j - 1))/ r1 - mu_xy * (full_dx(i, j) - full_dx(i - 1, j - 1));
            f_x = f_x - kxy*(r2 - lxy)* (full_x(i, j) - full_x(i - 1, j + 1))/ r2 - mu_xy * (full_dx(i, j) - full_dx(i - 1, j + 1));
            f_x = f_x + kxy*(r3 - lxy)* (-full_x(i, j) + full_x(i + 1, j - 1))/ r3 - mu_xy * (full_dx(i, j) - full_dx(i + 1, j - 1));
            f_x = f_x + kxy*(r4 - lxy)* (-full_x(i, j) + full_x(i + 1, j + 1))/ r4 - mu_xy * (full_dx(i, j) - full_dx(i + 1, j + 1));
            
            f_y = 0;
            % On y direction
            f_y = f_y - kxy*(r1 - lxy)* (full_y(i, j) - full_y(i - 1, j - 1))/ r1 - mu_xy * (full_dy(i, j) - full_dy(i - 1, j - 1));
            f_y = f_y - kxy*(r2 - lxy)* (full_y(i, j) - full_y(i - 1, j + 1))/ r2 - mu_xy * (full_dy(i, j) - full_dy(i - 1, j + 1));
            f_y = f_y + kxy*(r3 - lxy)* (-full_y(i, j) + full_y(i + 1, j - 1))/ r3 - mu_xy * (full_dy(i, j) - full_dy(i + 1, j - 1));
            f_y = f_y + kxy*(r4 - lxy)* (-full_y(i, j) + full_y(i + 1, j + 1))/ r4 - mu_xy * (full_dy(i, j) - full_dy(i + 1, j + 1));
            
            dstate_x(idx) = full_dx(i, j); 
            dstate_y(idx) = full_dy(i, j);
            ddstate_x(idx) = f_x/m; 
            ddstate_y(idx) = f_y/m;
            
            idx = idx + 1;
        end
    end
    
    f_x = [dstate_x, dstate_y, ddstate_x, ddstate_y]';
end