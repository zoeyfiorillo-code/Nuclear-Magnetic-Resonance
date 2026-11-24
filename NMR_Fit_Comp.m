%% Nuclear Magnetic Resonance - Program to compare exponential curve fitting to linear curve fitting
clc; clear all; close all;  help NMR_Fit_Comp      % Clear the memory and print header
%% read CSV file data scraped from oscilloscope
NMR_data = readtable('Data/NMR_Data.csv');
NMR_data.Properties.VariableNames = {'t_FID', 'V_FID'};
% extract columns and transpose into row vectors
t = NMR_data.t_FID';
V = NMR_data.V_FID';
%% Exponential curve fitting
% Exponential function
FID = @(b,t) b(1).*exp(-(t - 0.00047)./b(2))+b(3);
FID_b_0 = [1, 0.005, 0.24];   % [A0, T2*, C] initial guess
% Sigma for exponential
alpha = 0.01;
sigma_exp = alpha*ones(1,length(t));
% Error function and fitting for exponential
err_FID = @(b) sum(((V - FID(b,t))./sigma_exp).^2);
FID_b_fit = fminsearch(err_FID,FID_b_0);
A_0_fit = FID_b_fit(1);
T2_fit = FID_b_fit(2);
C=FID_b_fit(3);
V_fit = FID(FID_b_fit, t);
% Chi squared for exponential
res_exp = V - V_fit;
NM = length(t)-length(FID_b_fit);
chisqr_exp = sum((res_exp./sigma_exp).^2);
% Parameter uncertainty for exponential
J_FID = zeros(length(t), length(FID_b_fit));
db = 1e-6*abs(FID_b_fit) + 1e-8;
for j = 1:length(FID_b_fit)
    FID_bp = FID_b_fit;
    FID_bm = FID_b_fit;
    FID_bp(j) = FID_bp(j) + db(j);
    FID_bm(j) = FID_bm(j) - db(j);
    J_FID(:,j) = (FID(FID_bp,t) - FID(FID_bm,t))'/(2*db(j))./ sigma_exp';
end
Cov_FID = inv(J_FID'*J_FID);
SE_FID = sqrt(diag(Cov_FID));
A0_err = SE_FID(1);
T2_exp_err = SE_FID(2);
C_err = SE_FID(3);
% Print fit parameters for exponential
disp('Exponential fit parameters')
fprintf(' A0 = %.3f ± %.3f V\n', FID_b_fit(1), A0_err);
fprintf(' T2* = %.7f ± %.7f s\n', FID_b_fit(2), T2_exp_err);
fprintf(' C = %.3f ± %.3f V\n\n', FID_b_fit(3), C_err);
%% Plot exponential for FID T2*
figure(1); clf;
errorbar(t, V, sigma_exp,'.', 'MarkerSize', 10, 'Color', '#1171BE');
hold on
plot(t, V_fit, '-', 'LineWidth', 2, 'Color', '#DD5400');
xlabel('Time (sec)');
ylabel('Signal Height (Volts)');
legend('Original Data','Exponential Fit');
title(['Exponential Curve Fitting for T_2*' newline ...
    '\chi^2 = ', sprintf('%.0f    N-M = %.0f', chisqr_exp, NM)]);
% Equation and text box
Eq_1 = sprintf('$\\mathrm{Sig}(t)=%.3f\\,e^{-\\frac{t-0.00047}{%.4f}}%+.3f$', ...
                     FID_b_fit(1), FID_b_fit(2), FID_b_fit(3));
xpos = 0.0154; ypos = 1.2392;
text(xpos, ypos, Eq_1, ...
    'Interpreter','latex', ...
    'FontSize',14, ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'Margin',5);
grid on;
%% Linear curve fitting
lin_t = t - 0.00047;
lin_V = log((V-C)./2.38);
sigma = alpha./(V-C);
% Linear regression (Straight line) fit using linreg.m
% Author: Alejandro L. Garcia, Copyright 1999
% Numerical Methods for Physics, 2nd Ed. (Prentice Hall)
[a_fit_line sig_a_line yy_line chisqr_line] = linreg(lin_t,lin_V,sigma);     
% T2 extraction from fit
lin_T2 = -1/a_fit_line(2);
T2_lin_err = sig_a_line(2)/(a_fit_line(2)^2);
% Print fit parameters and T2* for linearized data
disp('Linearized data parameters')
fprintf(' a = %.2f ± %.2f \n', a_fit_line(2), sig_a_line(2));
fprintf(' b = %.4f ± %.4f \n', a_fit_line(1), sig_a_line(1));
fprintf(' T2* = %.7f ± %.7f s \n', lin_T2, T2_lin_err);
%% Plot linearized data for FID T2*
figure(2); clf;
errorbar(lin_t, lin_V, sigma,'.', 'MarkerSize', 10, 'Color', '#1171BE');
hold on
plot(lin_t, yy_line, '-', 'LineWidth', 2, 'Color', '#DD5400');
xlabel('t - 0.00047 (sec)');
ylabel('ln( Sig(t) / 2.38 )');
legend('Linearized Data', 'Linear Fit')
title(['Linear Curve Fitting for T_2*' newline ...
    '\chi^2 = ', sprintf('%.0f    N-M = %.0f', chisqr_line, NM)]);
% Equation and text box
Eq_2 = sprintf(['$\\mathrm{ln}\\hspace{0.1 mm}\\left(\\frac{\\mathrm{Sig}(t)}{2.38}\\right)=%.2f\\hspace{0.1 mm}\\left(\\hspace{0.1 mm}t-0.00047\\right) + \\left(%.4f\\right)$'], ...
                     a_fit_line(2), a_fit_line(1));
xpos = 0.0120; ypos = -1.0771;
text(xpos, ypos, Eq_2, ...
    'Interpreter','latex', ...
    'FontSize',14, ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'Margin',5);
grid on;
