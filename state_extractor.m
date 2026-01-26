function [x, y, dx, dy] = state_extractor(state)
    l = size(state, 1);
    n = l / 4;
    x = state(1: n, :);
    y = state(n + 1: 2*n, :);
    dx = state(2*n + 1: 3*n, :);
    dy = state(3*n + 1: 4*n, :);
end
