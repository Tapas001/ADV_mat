%% Loading data and taking gradient
M = xlsread('AFNew_vel_data.csv');
y = M(:,1);
U = M(:,2);
u_grad = gradient(U,y);
%% this one is for equation 3.1
x = 1./y;
y = u_grad;
m = sum(u_grad.*x)/sum(x.^2); %when the y-intercept is zero this is the formula for slope
kappa = 0.41; % von-kaarman constant
u_starfirst = kappa*m;

%% this one is for equation 3.2
z1 = M(:,1);
U = M(:,2);
z = z1+0.0133; % here change the value of 0.0133 for different data set to match u_starfirst and u_starsecond
%
y1 = log(z);

% Calculate the linear regression coefficients and r2 value
[a0, a1, r2] = calculate_linear_regression_coef(y1,U );

u_starsecond = a1 * kappa; % the value of shear velocity from second equation

    % Number of observations/points
    n = numel(x);

    sx = sum(x);
    sy = sum(y);

    sx2 = dot(x, x);
    sxy = dot(x, y);
    sy2 = dot(y, y);

    % Calculating regression coefficients
    a1 = (n * sxy - sx * sy) / (n * sx2 - sx^2);
    a0 = sy / n - a1 * sx / n;

    r2 = ((n * sxy - sx * sy) / sqrt(n * sx2 - sx^2) / sqrt(n * sy2 - sy^2))^2;
end
