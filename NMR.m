%% NMR - Program to fit NMR Data
clc; clear all; close all; help NMR      % Clear the memory and print header
%% T2* Effective Transverse Relaxation Time via FID
% read CSV file data scraped from oscilloscope
NMR_data = readtable('NMR_Data.csv');
NMR_data.Properties.VariableNames = {'t_FID', 'V_FID'};
% extract columns and transpose into row vectors
t_FID = NMR_data.t_FID';
V_FID = NMR_data.V_FID';
f_FID = @(b_FID,t_FID) b_FID(1).*exp(-(t_FID - 0.00048)./b_FID(2))+b_FID(3);
b_FID_0 = [1, 0.005, 0.24];     % [A0, T2*, C] initial guess
err_FID = @(b_FID) sum((V_FID - f_FID(b_FID,t_FID)).^2);
b_FID_fit = fminsearch(err_FID,b_FID_0);
A_0_fit = b_FID_fit(1);
T2_FID_fit = b_FID_fit(2);
C_FID_fit = b_FID_fit(3);
t_FID_fit = linspace(min(t_FID),max(t_FID),300);
V_FID_fit = f_FID(b_FID_fit,t_FID_fit);
V_FID_fit_data = f_FID(b_FID_fit, t_FID);
SS_res_FID = sum((V_FID - V_FID_fit_data).^2);
SS_tot_FID = sum((V_FID - mean(V_FID)).^2);
R2_FID = 1 - SS_res_FID/SS_tot_FID;
%% T2* parameter uncertainty for FID
res_FID = V_FID - V_FID_fit_data;
n_FID = length(V_FID);      
p_FID = length(b_FID_fit);
% residual variance
s2_FID = sum(res_FID.^2) / (n_FID - p_FID);
% finite-difference Jacobian
J_FID = zeros(n_FID, p_FID);
db_FID = 1e-6 * abs(b_FID_fit) + 1e-8;
for j = 1:p_FID
    b_FID_p = b_FID_fit;
    b_FID_m = b_FID_fit;
    b_FID_p(j) = b_FID_p(j) + db_FID(j);
    b_FID_m(j) = b_FID_m(j) - db_FID(j);
    J_FID(:,j) = (f_FID(b_FID_p, t_FID) - f_FID(b_FID_m, t_FID))' / (2*db_FID(j));
end
con_FID = cond(J_FID' * J_FID);
% covariance matrix and standard errors
Cov_FID = s2_FID * pinv(J_FID' * J_FID);
SE_FID = sqrt(diag(Cov_FID));
A0_FID_err = SE_FID(1);
T2_FID_err = SE_FID(2);
C_FID_err = SE_FID(3);
fprintf('A_0 = %.4f ± %.4f V\n', b_FID_fit(1), A0_FID_err);
fprintf('T2* = %.6f ± %.6f s\n', b_FID_fit(2), T2_FID_err);
fprintf('C = %.4f ± %.4f V\n', b_FID_fit(3), C_FID_err);
fprintf('Condition number = %f\n', con_FID);
%% Plot for T2* for FID
figure(1); clf;
plot(t_FID, V_FID, '.', 'MarkerSize', 10, 'Color', '#1171BE');
hold on
plot(t_FID_fit, V_FID_fit, '-', 'LineWidth', 2, 'Color', '#DD5400');
xlabel('t (s)');
ylabel('Signal Height (V)');
legend('Data','Fit');
title('T_2* Decay Time Via FID');
Eq_T2_FID = sprintf(['$\\mathrm{Sig}(t)=%.3f\\,e^{-\\frac{t-0.00048}{%.4f}}%+.3f$' newline ...
                     '$R^{2}=%.3f$'], ...
                     b_FID_fit(1), b_FID_fit(2), b_FID_fit(3), R2_FID);
xpos = 0.0153; ypos = 1.2083;
text(xpos, ypos, Eq_T2_FID, ...
    'Interpreter','latex', ...
    'FontSize',14, ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'Margin',5);
grid on;
%% T2 Spin-Spin Relaxation Time via echo heights
t2 = [0.0150 0.0200 0.0250 0.0300 0.0350 0.0400 0.0450 0.0500 0.0550 0.0600 0.0650 0.0700 0.0750 0.0800 0.0850 0.0900 0.0950 0.1000];
V2 = [1.210 1.130 1.060 1.000 0.928 0.872 0.792 0.744 0.672 0.624 0.560 0.496 0.448 0.392 0.344 0.296 0.256 0.208];
f2 = @(b2,t2) b2(1).*exp(-2*t2./b2(2));
b2_0 = [1, 0.150];                        % [M0, T2] initial guess
err2 = @(b2) sum((V2 - f2(b2,t2)).^2);
b2_fit = fminsearch(err2,b2_0);
M2_0_fit = b2_fit(1);
T2_fit = b2_fit(2);
t2fit = linspace(min(t2),max(t2),300);
V2fit = f2(b2_fit,t2fit);
V2fit_data = f2(b2_fit, t2);
SS_res2 = sum((V2 - V2fit_data).^2);
SS_tot2 = sum((V2 - mean(V2)).^2);
R2_T2 = 1 - SS_res2/SS_tot2;
%% T2 parameter uncertainty for echo
res2 = V2 - V2fit_data;
n2 = length(V2);      
p2 = length(b2_fit);
% residual variance
s2_2 = sum(res2.^2) / (n2 - p2);
% finite-difference Jacobian
J2 = zeros(n2, p2);
db2 = 1e-6 * abs(b2_fit) + 1e-8;
for j = 1:p2
    b2p = b2_fit;
    b2m = b2_fit;
    b2p(j) = b2p(j) + db2(j);
    b2m(j) = b2m(j) - db2(j);
    J2(:,j) = (f2(b2p, t2) - f2(b2m, t2))' / (2*db2(j));
end
con_T2 = cond(J2' * J2);
% covariance matrix and standard errors
Cov2 = s2_2 * inv(J2' * J2);
SE2 = sqrt(diag(Cov2));
M0_T2_err = SE2(1);
T2_err = SE2(2);
fprintf('\nM_0 = %.4f ± %.4f V\n', b2_fit(1), M0_T2_err);
fprintf('T2 = %.4f ± %.4f s\n', b2_fit(2), T2_err);
fprintf('Condition number = %f\n', con_T2);
%% Plot for T2 for echo
figure(2); clf;
plot(t2, V2, '.', 'MarkerSize', 15, 'Color', '#1171BE');
hold on
plot(t2fit, V2fit, '-','LineWidth', 2, 'Color', '#DD5400');
xlabel('\tau (s)');
ylabel('Signal Height (V)');
legend('Data','Fit');
title('T_2 Spin-Spin Relaxation Time Via Spin-Echo');
Eq_T2 = sprintf(['$M_{xy}(2\\tau)=%.3f\\,e^{-\\,\\frac{2\\tau}{%.4f}}$' newline ...
                 '$R^{2}=%.3f$'], ...
                  b2_fit(1), b2_fit(2), R2_T2);
xpos = 0.0843; ypos = 1.2120;
text(xpos, ypos, Eq_T2, ...
    'Interpreter','latex', ...
    'FontSize',14, ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'Margin',5);
grid on;
%% T1 Spin-Lattice Relaxation Time
t1 = [0.0101 0.0200 0.0300 0.0400 0.0500 0.0600 0.0700 0.0800 0.0900 0.1000 0.11100 0.1200 0.1300 0.1400 0.1500 0.1600];
V1 = [1.330 1.230 1.150 1.050 0.952 0.856 0.768 0.664 0.57 0.472 0.384 0.288 0.192 0.112 0.056 0.000];
f1 = @(b1,t1) abs(b1(1).*(1 - 2*exp(-t1./b1(2))));  
b1_0 = [1, 0.153];                                  % [M0, T1] initial guess
err1 = @(b1) sum((V1 - f1(b1,t1)).^2);
b1_fit = fminsearch(err1,b1_0);
M1_0_fit = b1_fit(1);
T1_fit = b1_fit(2);
t1fit = linspace(min(t1),max(t1),300);
V1fit = f1(b1_fit,t1fit);
V1fit_data = f1(b1_fit, t1);
SS_res1 = sum((V1 - V1fit_data).^2);
SS_tot1 = sum((V1 - mean(V1)).^2);
R2_T1 = 1 - SS_res1/SS_tot1;
%% T1 parameter uncertainty
res1 = V1 - V1fit_data;
n1 = length(V1);      % number of data points
p1 = length(b1_fit);  % number of fitted parameters
% residual variance
s2_1 = sum(res1.^2) / (n1 - p1);
% finite-difference Jacobian
J1 = zeros(n1, p1);
db1 = 1e-6 * abs(b1_fit) + 1e-8;
for j = 1:p1
    b1p = b1_fit;
    b1m = b1_fit;
    b1p(j) = b1p(j) + db1(j);
    b1m(j) = b1m(j) - db1(j);
    J1(:,j) = (f1(b1p, t1) - f1(b1m, t1))' / (2*db1(j));
end
con_T1 = cond(J1' * J1);
% covariance matrix and standard errors
Cov1 = s2_1 * inv(J1' * J1);
SE1 = sqrt(diag(Cov1));
M0_T1_err = SE1(1);
T1_err = SE1(2);
fprintf('\nM_0 = %.4f ± %.4f V\n', b1_fit(1), M0_T1_err);
fprintf('T1 = %.4f ± %.4f s\n', b1_fit(2), T1_err);
fprintf('Condition number = %f\n', con_T1);
%% Plot for T1
figure(3); clf;
plot(t1, V1, '.', 'MarkerSize', 15, 'Color', '#1171BE');
hold on
plot(t1fit, V1fit, '-', 'LineWidth', 2, 'Color', '#DD5400');
xlabel('t (s)');
ylabel('Signal Height (V)');
legend('Data','Fit');
title('T_1 Spin-Lattice Relaxation Time');
Eq_T1 = sprintf(['$\\left|M_z(t)\\right|=%.3f\\left|1-2e^{-\\,\\frac{t}{%.4f}}\\right|$' newline ...
                 '$R^{2}=%.3f$'], ...
                  b1_fit(1), b1_fit(2), R2_T1);
xpos = 0.1269; ypos = 1.1866;
text(xpos, ypos, Eq_T1, ...
    'Interpreter','latex', ...
    'FontSize',14, ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'Margin',5)
grid on;
%% To get original graph colors
RGB = orderedcolors("gem");
H = rgb2hex(RGB);