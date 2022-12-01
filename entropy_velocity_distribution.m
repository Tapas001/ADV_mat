clc;
data = xlsread('data.xlsx');
%% default plot settings
set(groot,'defaultLineLineWidth',2.5) 
set(0,'DefaultaxesLineWidth', 1.5) 
set(0,'DefaultaxesFontSize', 14) 
set(0,'DefaultaxesFontWeight', 'bold') 
set(0,'DefaultTextFontSize', 14) 
set(0,'DefaultaxesFontName', 'latex') 
set(0,'DefaultlegendFontName', 'latex')
set(0,'defaultAxesXGrid','off') 
set(0,'defaultAxesYGrid','off') 
set(0,'DefaultTextInterpreter','latex') 


%read y
z = data(1:end,1);

%D = max(z); %depth
D = max(z)+0.02; %depth % 5cn is added for ADV limitation

%read velocity
u = data(1:end,2);

[umax, index_of_umax] = max(u);
delta = z(index_of_umax);
p(1) = 5.0;

%test = umax/p(1) * log( 1 + (exp(p(1))-1)*(z/delta).*exp(1-z/delta));

% Define your function
fh = @(x,p,umax,delta) umax/p(1) * log( 1 + (exp(p(1))-1)*(x/delta).*exp(1-x/delta));

errfh = @(p,x,y,umax,delta) sum((y(:)-fh(x(:),p,umax,delta)).^2);

% an initial guess of the M
p0 = [5.0];
    
 % search for solution
P1 = fminsearch(errfh,p0,[],z,u,umax,delta);

comp1 = fh(z,P1,umax,delta);

%compute R2 value
SStot = sum((u - mean(comp1)).^2);      % Compute R-squared
SSres = errfh(P1,z,u,umax,delta);
Rsq   = 1 - (SSres/SStot);

display('M=');
display(P1(1));
display('Umax');
display(umax);
display('h');
display(D-delta);
display('depth');
display(D);
display('R2 value=');
display(Rsq);

plot(u/umax,z/D,'ro');
hold on
plot(comp1/umax,z/D);
hold off
axis([0 1.2 0 1.2]);
% axis([0 0.6 0 0.2]);
legend('Measured','Computed');
xlabel('u/umax') 
ylabel('y/D')

