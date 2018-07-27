function [ BC ] = bhattacharyya( Vector1, Vector2 )
%bhattacharya Calculates Bhattacharyya distance between time series of two
%pixels and Bhattacharyya coefficient
%   Code in R taken from here:
%   http://stats.stackexchange.com/questions/78849/measure-for-separability/78855#78855

% convert vectors to matrices in case they are not
Matrix1 = vec2mat(Vector1,1);
Matrix2 = vec2mat(Vector2,1);

% define means
mean_Matrix1 = mean(Matrix1);
mean_Matrix2 = mean(Matrix2);

% define difference of means
mean_difference = mean_Matrix1 - mean_Matrix2;

% define covariances for supplied matrices
cv_Matrix1 = cov(Matrix1);
cv_Matrix2 = cov(Matrix2);

% define the halfsum of cv's as "p"
p = (cv_Matrix1 + cv_Matrix2)./2;

%% distance and coefficient calculation 

% calculate the Bhattacharyya index
bh_distance = 0.125 * (mean_difference) * p^(-1) * mean_difference + 0.5...
    * log(det(p)/sqrt(det(cv_Matrix1) * det(cv_Matrix2)));

% actually want Bhattacharyya coefficient, because that is 
% limited between 0 and 1
% the distance is the negative natural log of the coefficient
BC = 1./exp(bh_distance);

end