% Load the data from file
data = xlsread('velocity_data1.xlsx'); % put your depth wise data for caiculating the velocity profile
z = data(:, 1);
U = data(:, 2);

% Calculate the logarithm of z
x = log(z);
y = U;

% Calculate the linear regression coefficients and r2 value
[a0, a1, r2] = calculate_linear_regression_coef(x, y);

kappa = 0.41;
u_star = a1 * kappa;
fprintf('Fitted shear velocity is %.3f m/s.\n', u_star);

ks = exp(a0 * kappa / (-u_star)) * 30;
fprintf('Fitted roughness is %.4f m.\n', ks);

U_pre = (u_star / kappa) * log(z / ks * 30);

% Plot the data and fitted curve
scatter(U, z, 'r', 'x', 'LineWidth', 1.5, 'DisplayName', 'experimental data');
hold on;
plot(U_pre, z, 'LineWidth', 1.5, 'DisplayName', 'fitted curve');
hold off;

xlabel('U (m/s)', 'FontSize', 16);
ylabel('z (m)', 'FontSize', 16);

% Set the font size of ticks
set(gca, 'FontSize', 12);

% Add legend
legend('Location', 'northwest', 'FontSize', 14, 'AutoUpdate', 'off');

% Function to calculate linear regression coefficients and r2 value
function [a0, a1, r2] = calculate_linear_regression_coef(x, y)
    if numel(x) ~= numel(y)
        error('The two vectors x and y are not of the same length.');
    end

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
