% X Ϊ n*p ����, ���� n Ϊ�۲���, p Ϊά��
% y Ϊ n*1 ����, ȡֵΪ {0, 1}
clear; clc
data=csvread('data.csv');
X=data(:,2);
y=data(:,3);
n=length(y);
N=59999;%��������
lt=2; %��������
beta=[zeros(N+1,2) ones(N+1,2)];
average=zeros(N+1,2*lt);
X_ext = [ones(n, 1), X];

for i=1:lt
    for t=2:N+1
        logratio=X_ext*beta(t-1,lt*(i-1)+1:lt*(i-1)+2)';
        u=rand(n,1).*exp(y.*logratio)./(1+exp(logratio));
        u_y1=u(1:14);
        u_y0=u(15:54);
        L0=max(log(u_y1./(1-u_y1))-beta(t-1,lt*(i-1)+2).*X(1:14));
        U0=min(log((1-u_y0)./u_y0)-beta(t-1,lt*(i-1)+2).*X(15:54));
        beta(t,lt*(i-1)+1)=100*trandn(L0/100,U0/100);
        L1=max((log(u_y1./(1-u_y1))-beta(t,lt*(i-1)+1))./X(1:14));
        U1=min((log((1-u_y0)./u_y0)-beta(t,lt*(i-1)+1))./X(15:54));
        beta(t,lt*(i-1)+2)=100*trandn(L1/100,U1/100);
        average(t,lt*(i-1)+1:lt*(i-1)+2)=(average(t-1,lt*(i-1)+1:lt*(i-1)+2)*(t-1)+beta(t,lt*(i-1)+1:lt*(i-1)+2))/t;
    end
end
%%
%��Ե�һ����
fprintf('beta0�ĺ�����������Ϊ%.4f \n',(average(N+1,1)*(N+1)-average(25000,1)*25000)/(N+1-25000));
fprintf('beta1�ĺ�����������Ϊ%.4f \n',(average(N+1,2)*(N+1)-average(25000,2)*25000)/(N+1-25000));
%HPD�������� monte carlo��
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
fprintf('beta0��HPD��������Ϊ[%.4f,%.4f] \n',beta0_l,beta0_u);
fprintf('beta1��HPD��������Ϊ[%.4f,%.4f] \n',beta1_l,beta1_u);
%����·��
figure(1);
plot(beta(:,1))
ylabel('beta0');title('����·��');
figure(2);
plot(beta(:,2))
ylabel('beta1');title('����·��');
%������ֵ
figure(3);
plot(average(:,1))
ylabel('beta0');title('������ֵ');
figure(4);
plot(average(:,2))
ylabel('beta1');title('������ֵ');
%ɢ��
figure(5);
scatter(beta(25001:N+1,1),beta(25001:N+1,2),'.')
xlabel('beta0');ylabel('beta1');title('����ֲ�ɢ��ͼ');
%ֱ��ͼ
figure(6);
subplot(1,2,1);
histogram(beta(25001:N+1,1));
xlabel('beta0');title('��������ֱ��ͼ');
subplot(1,2,2);
histogram(beta(25001:N+1,2));
xlabel('beta1');title('��������ֱ��ͼ');
%%
%GR��
beta0_fai1=mean(beta(:,1));
beta0_fai2=mean(beta(:,3));
beta0_fai=mean(mean(beta(:,[1,3])));
beta0_Bn=(beta0_fai1-beta0_fai)^2+(beta0_fai2-beta0_fai)^2;
beta0_s1=var(beta(:,1));
beta0_s2=var(beta(:,3));
beta0_W=(beta0_s1+beta0_s2)/lt;
beta0_V=(n-1)*beta0_W/n+beta0_Bn+beta0_Bn/lt;
beta0_R=beta0_V/beta0_W;
fprintf('beta0��R��ֵֵΪ%.4f \n',beta0_R);

beta1_fai1=mean(beta(:,2));
beta1_fai2=mean(beta(:,4));
beta1_fai=mean(mean(beta(:,[2,4])));
beta1_Bn=(beta1_fai1-beta1_fai)^2+(beta1_fai2-beta1_fai)^2;
beta1_s1=var(beta(:,2));
beta1_s2=var(beta(:,4));
beta1_W=(beta1_s1+beta1_s2)/lt;
beta1_V=(n-1)*beta1_W/n+beta1_Bn+beta1_Bn/lt;
beta1_R=beta1_V/beta1_W;
fprintf('beta1��R��ֵֵΪ%.4f \n',beta1_R);