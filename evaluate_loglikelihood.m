function loglikelihood = evaluate_loglikelihood(X, y, beta)
% x Ϊ n*p ����, ���� n Ϊ�۲���, p Ϊά��
% y Ϊ n*1 ����, ȡֵΪ {0, 1}
% beta Ϊ p*1 ϵ������

n = size(X, 1);
X_ext = [ones(n, 1), X];
logratio = X_ext * beta;
zeros_x = zeros(n, 1);
loglikelihood = - sum( y.*logsumexp([zeros_x, -logratio], 2) + (1 - y).*logsumexp([zeros_x, logratio], 2) );