function MSE = wmse(residual, weight)
% Weighted mean squared error. If weight is interpreted as probability this
% corresponds the expected value of the square of residual.
arguments
    residual double
    weight double {mustBePositive}
end

MSE = sum(weight.*residual.^2, 'all')/sum(weight, 'all');
end