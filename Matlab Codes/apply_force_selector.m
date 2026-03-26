function [B_mod, u_mod] = apply_force_selector(B_full, u_full, force_selector)
% APPLY_FORCE_SELECTOR
% Selects columns of B_full corresponding to forces applied to selected masses.
%
% Inputs:
%   B_full          : symbolic model (size nz x 2*N)
%   u_full          : symbolic inout (size nz x 2*N)
%   force_selector  : logical or 0/1 vector of length N
%
% Output:
%   B_mod          : modified B with only active forces (2 columns per selected mass)
%   u_mod          : modified u with only active forces (2 columns per selected mass)
%
% Example:
%   force_selector = [1 0 1]; % apply forces on mass 1 and 3
%   [B_mod, u_mod] = apply_force_selector(model, force_selector);

N = length(force_selector);
cols_to_keep = [];

shift = size(B_full, 2) / 2;

for i = 1:N
    if force_selector(i)
        idx = i;
        cols_to_keep = [cols_to_keep, idx, idx+shift]; % x and y columns
    end
end

B_mod = B_full(:, cols_to_keep);
u_mod = u_full(cols_to_keep);
end