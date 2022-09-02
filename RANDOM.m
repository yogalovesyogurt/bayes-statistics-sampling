% X 为 n*p 矩阵, 其中 n 为观测数, p 为维数
% y 为 n*1 向量, 取值为 {0, 1}
clear; clc
data=csvread('data.csv');
X=data(:,2);
y=data(:,3);
n=length(y);
N=59999;%迭代次数
lt=2; %链的条数
beta=[zeros(N+1,2) ones(N+1,2)];
average=zeros(N+1,2*lt);
count=zeros(lt,1);
for i=1:lt
    for t=2:N+1
        beta(t,lt*(i-1)+1)=normrnd(beta(t-1,lt*(i-1)+1),0.1);
        beta(t,lt*(i-1)+2)=normrnd(beta(t-1,lt*(i-1)+2),0.1);
        u=rand;
        nume=evaluate_loglikelihood(X, y, beta(t,lt*(i-1)+1:lt*(i-1)+2)')+log(normpdf(beta(t,lt*(i-1)+1),0,100)*normpdf(beta(t,lt*(i-1)+2),0,100));
        deno=evaluate_loglikelihood(X, y, beta(t-1,lt*(i-1)+1:lt*(i-1)+2)')+log(normpdf(beta(t-1,lt*(i-1)+1),0,100)*normpdf(beta(t-1,lt*(i-1)+2),0,100));
        alpha=min(1,exp(nume-deno));
        if u>alpha
            beta(t,lt*(i-1)+1:lt*(i-1)+2)=beta(t-1,lt*(i-1)+1:lt*(i-1)+2);
        else
            count(i)=count(i)+1;
        end
        average(t,lt*(i-1)+1:lt*(i-1)+2)=(average(t-1,lt*(i-1)+1:lt*(i-1)+2)*(t-1)+beta(t,lt*(i-1)+1:lt*(i-1)+2))/t;
    end
end
%%
%针对第一条链
accept=count(1)/N; %接受率
fprintf('接受概率为%.4f \n',accept);
fprintf('beta0的后验期望估计为%.4f \n',(average(N+1,1)*(N+1)-average(25000,1)*25000)/(N+1-25000));
fprintf('beta1的后验期望估计为%.4f \n',(average(N+1,2)*(N+1)-average(25000,2)*25000)/(N+1-25000));
%HPD置信区间 monte carlo法
p=0.9;
sorted_beta0=sort(beta(25001:N+1,1));
sorted_beta1=sort(beta(25001:N+1,2));
n_in_ci=ceil(p*(N+1-25000));
n_ci=N+1-25000-n_in_ci;
ci_width=sorted_beta0((1:n_ci)+n_in_ci)-sorted_beta0(1:n_ci);
[~,idx_1]=min(ci_width);
beta0_l=sorted_beta0(idx_1);
beta0_u=sorted_beta0(idx_1+n_in_ci);
ci_width1=sorted_beta1((1:n_ci)+n_in_ci)-sorted_beta1(1:n_ci);
[~,idx_2]=min(ci_width1);
beta1_l=sorted_beta1(idx_2);
beta1_u=sorted_beta1(idx_2+n_in_ci);
fprintf('beta0的HPD可信区间为[%.4f,%.4f] \n',beta0_l,beta0_u);
fprintf('beta1的后验期望估计为[%.4f,%.4f] \n',beta1_l,beta1_u);
%样本路径
figure(1);
plot(beta(:,1))
ylabel('beta0');title('样本路径');
figure(2);
plot(beta(:,2))
ylabel('beta1');title('样本路径');
%遍历均值
figure(3);
plot(average(:,1))
ylabel('beta0');title('遍历均值');
figure(4);
plot(average(:,2))
ylabel('beta1');title('遍历均值');
%散点
figure(5);
scatter(beta(25001:N+1,1),beta(25001:N+1,2),'.')
xlabel('beta0');ylabel('beta1');title('后验分布散点图');
%直方图
figure(6);
subplot(1,2,1);
histogram(beta(25001:N+1,1));
xlabel('beta0');title('后验样本直方图');
subplot(1,2,2);
histogram(beta(25001:N+1,2));
xlabel('beta1');title('后验样本直方图');
%%
%GR法
beta0_fai1=mean(beta(:,1));
beta0_fai2=mean(beta(:,3));
beta0_fai=mean(mean(beta(:,[1,3])));
beta0_Bn=(beta0_fai1-beta0_fai)^2+(beta0_fai2-beta0_fai)^2;
beta0_s1=var(beta(:,1));
beta0_s2=var(beta(:,3));
beta0_W=(beta0_s1+beta0_s2)/lt;
beta0_V=(n-1)*beta0_W/n+beta0_Bn+beta0_Bn/lt;
beta0_R=beta0_V/beta0_W;
fprintf('beta0的R估值值为%.4f \n',beta0_R);

beta1_fai1=mean(beta(:,2));
beta1_fai2=mean(beta(:,4));
beta1_fai=mean(mean(beta(:,[2,4])));
beta1_Bn=(beta1_fai1-beta1_fai)^2+(beta1_fai2-beta1_fai)^2;
beta1_s1=var(beta(:,2));
beta1_s2=var(beta(:,4));
beta1_W=(beta1_s1+beta1_s2)/lt;
beta1_V=(n-1)*beta1_W/n+beta1_Bn+beta1_Bn/lt;
beta1_R=beta1_V/beta1_W;
fprintf('beta1的R估值值为%.4f \n',beta1_R);