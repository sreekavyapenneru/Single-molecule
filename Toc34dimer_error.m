% Given data
d = [49.3, 35.4, 24.3];  % Dimer values
m = [50.7	64.57512433	75.70405554];  % Monomer values

% Compute f values
f = d ./ (d + m);

% Mean of f
mean_f = mean(f);

% Standard deviations of d and m
sigma_d = std(d, 1);  % Use 1 for population standard deviation
sigma_m = std(m, 1);

% Compute covariance between d and m
cov_dm = cov(d, m);
sigma_dm = cov_dm(1, 2);  % Extract covariance value

% Compute partial derivatives
denom = (d + m).^2;
df_dd = m ./ denom;
df_dm = -d ./ denom;

% Propagation of uncertainty formula
sigma_f_squared = (df_dd.^2 * sigma_d^2) + (df_dm.^2 * sigma_m^2) + 2 * (df_dd .* df_dm * sigma_dm);
sigma_f = sqrt(sum(sigma_f_squared) / length(d));  % Standard error of f

% Display results
fprintf('Mean f: %.4f\n', mean_f);
fprintf('Standard Error of f: %.4f (%.2f%%)\n', sigma_f, sigma_f * 100);
