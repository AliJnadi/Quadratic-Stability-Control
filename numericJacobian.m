function J = numericJacobian(fhandle, x)
%NUMERICJACOBIAN
% Central finite-difference approximation of Jacobian

    n = numel(x);
    J = zeros(n,n);

    % Adaptive step size
    h = sqrt(eps) * max(1, norm(x));

    fx = fhandle(x);
    if numel(fx) ~= n
        error('numericJacobian:DimensionMismatch', ...
              'Function output dimension does not match input.');
    end

    for i = 1:n
        dx = zeros(n,1);
        dx(i) = h;

        f_plus  = fhandle(x + dx);
        f_minus = fhandle(x - dx);

        J(:,i) = (f_plus - f_minus) / (2*h);
    end
end

% function J = numericJacobian(fhandle, x)
%     n = length(x);
%     J = zeros(n,n);
%     h = 1e-6;
% 
%     for i = 1:n
%         dx = zeros(n,1);
%         dx(i) = h;
% 
%         f_plus  = fhandle(x + dx);
%         f_minus = fhandle(x - dx);
% 
%         J(:,i) = (f_plus - f_minus) / (2*h);
%     end
% end