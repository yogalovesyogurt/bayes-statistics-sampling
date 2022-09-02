function loglikelihood = evaluate_loglikelihood(X, y, beta)
% x 为 n*p 矩阵, 其中 n 为观测数, p 为维数
% y 为 n*1 向量, 取值为 {0, 1}
% beta 为 p*1 系数向量

n = size(X, 1);
X_ext = [ones(n, 1), X];
logratio = X_ext * beta;
zeros_x = zeros(n, 1);
loglikelihood = - sum( y.*logsumexp([zeros_x, -logratio], 2) + (1 - y).*logsumexp([zeros_x, logratio], 2) );