function CCCs = concordanceSimilarity(x,y)
% x: column vector of samples from random variable X
% y: column vector of samples from random variable Y, of same length as x
% returns: concordance correlation coeffcient (CCC) between x and y
x = x(:); % flaten x if it is a matrix
y = y(:); % flaten y if it is a matrix
covariance = cov(x, y);
covariance = covariance(2,1);

CCCs = 2*covariance / (var(x) + var(y) + (mean(x) - mean(y)).^2);


end