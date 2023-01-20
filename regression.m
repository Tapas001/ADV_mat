%% Getting data


m = xlsread('t_regression2.xlsx');
k = m(:,1);
k1=m(:,2);
k2=m(:,3);
k3=m(:,4);
k4=m(:,5);
k5=m(:,6);
% Taking the logarith of the data will transfer the nonlinear system into
% multivariate linear system
y =log(m(:,1));
t1=log(m(:,2));
t2=log(m(:,3));
t3=log(m(:,4));
t4=log(m(:,5));
t5=log(m(:,6));
%% linear least square regression
n = length(y);
x= [ones(n,1) , t1,t2,t3,t4,t5];
p = y;
phi = inv(x'*x)*x'*y;
% Parameters
anot = exp(phi(1,1));
a1=phi(2,1);
a2 = phi(3,1);
a3 = phi(4,1);
a4 = phi(5,1);
a5 = phi(6,1);
% anot=exp(-1.588);
% a1=0.01;
% a2=-0.1367;
% a3=0.2584;
% a4=0.5129;
% a5=-0.934;

%% Statistical parameters
y_predicted = anot*[((k1).^a1).*((k2).^a2).*((k3).^a3).*((k4).^a4).*((k5).^a5)];
% RMSE
RMSE = sqrt((sum((y_predicted - k).^2)/(n)))
% MSE
MSE = (sum((y_predicted - k).^2)/(n))
% MAE 
MAE = sum(abs((k - y_predicted)))/(n)
% MAPE
MAPE = (100/n)*sum(abs(k - y_predicted)./(abs(k)));
% R_SQUARE
R_SQUARE = 1-(sum((y_predicted - k).^2))/(sum((k - mean(k)).^2))


%% Plotting
plot(y_predicted,k,"k*")
hold on
plot([0 1.5], [0 2])
hold on
plot([0 2], [0 1.5])
hold on
plot([0 2], [0 1.8])
hold on
plot([0 1.8], [0 2])
hold on
hold on
plot([0 2], [0 2])
%% Sensitivity Analysis
%% UNCERTAINITY ANALYSIS

