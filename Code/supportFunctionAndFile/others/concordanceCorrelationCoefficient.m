function ccc = concordanceCorrelationCoefficient(x,y)
% x: column vector of samples from random variable X
% y: column vector of samples from random variable Y, of same length as x
% returns: concordance correlation coeffcient (CCC) between x and y

covariance = cov(x, y);
covariance = covariance(2,1);

ccc = 2*covariance / (var(x) + var(y) + (mean(x) - mean(y)).^2);
end