%**************************************************************************
% This script runs examples of the Matlab wrapper functions to the main
% Madym C++ toolkit. It can be used to test the C++ tools are working as
% they should on your computer. However please refrain from modifying this
% script to the specifics of your computer (by all means make changes, but
% try to make them useful for everyone). If you need to change anything
% specific to your own set-up, please take a copy of the script, or
% generate the appropriate variables elswhere to be called from here
%--------------------------------------------------------------------------

%Set ground truth values for T1, S0 and 2CXM model parameters
n_vox = 24;
n_times = 100;

T1 = 1000 + randn(n_vox,1)*200;
S0 = 2000 + randn(n_vox,1)*500;
F_p = 1.0 + randn(n_vox,1)*0.1;
PS = 0.15 + randn(n_vox,1)*0.015;
v_e = 0.2 + randn(n_vox,1)*0.02;
v_p = 0.1 + randn(n_vox,1)*0.01;
    
TR = 3.5;
FAs = [2 10 18];
injection_img = 8;
dyn_t = linspace(0, 5, n_times);

dyn_FA = 20;
r1Const = 3.4;

%Generate VFA signals and add some noise
vfa_signals = signal_from_T1(T1, S0, FAs, TR);
vfa_signals = vfa_signals + randn(n_vox,3)*1;

%Generate AIF and modelled dynamic concentration time-series
Ca_t = population_aif(dyn_t, injection_img);
C_t = two_cx_model(F_p, PS, v_e, v_p,  0, Ca_t, dyn_t);

%Add some noise and rescale so baseline mean is 0
C_tn = C_t + randn(n_vox,n_times)*0.01;
try
    C_tn = C_tn - mean(C_tn(:,1:injection_img),2);
catch
    C_tn = bsxfun(@minus, C_tn, mean(C_tn(:,1:injection_img),2));
end

%Convert the concentrations to signals
S_tn = concentration_to_signal(...
    C_tn, dyn_FA, TR, T1, S0, r1Const, false);
%%
% 1) Fit baseline T1 from the VFA signals

%First run this in data mode using calculate_T1_lite:    
[T1_fit, S0_fit] = run_madym_T1(...
    'FAs', FAs,... 
    'signals', vfa_signals,...
    'TR', TR,... 
    'method', 'VFA',...
    'dummy_run', 0);
vfa_signals_fit = signal_from_T1(T1_fit, S0_fit, FAs, TR);

min_t = min(min(T1), min(T1_fit));
max_t = max(max(T1), max(T1_fit));
min_s = min(min(S0), min(S0_fit));
max_s = max(max(S0), max(S0_fit));
%%
figure('Name', 'Baseline T1');
subplot(1,2,1);
plot([min_t max_t], [min_t max_t], 'k--');
hold on;
plot(T1, T1_fit, 'rx');
title('T_1');
xlabel('Actual');
ylabel('Predicted');
axis equal;
axis([min_t max_t min_t max_t]);

subplot(1,2,2);
plot([min_s max_s], [min_s max_s], 'k--');
hold on;
plot(S0, S0_fit, 'rx');
title('S_0');
xlabel('Actual');
ylabel('Predicted');
axis equal;
axis([min_s max_s min_s max_s]);
%%
% 2) Fit the concentration curves using madym_lite
% a) Using the 2CXM
[model_params_C_2CXM, model_fit_C_2CXM, ~,~,CmC_t_2CXM] = run_madym_lite(...   
    '2CXM', C_tn, 'dyn_times', dyn_t);
% b) Using the 2CXM
[model_params_C_ETM, model_fit_C_ETM, ~,~,CmC_t_ETM] = run_madym_lite(...
    'ETM', C_tn, 'dyn_times', dyn_t);

display_dce_model_fits(C_tn, cat(3, CmC_t_2CXM, CmC_t_ETM), ...
    'plot_rows', 4, 'plot_cols', 6);
%%
% Look at parameter estimates for 2CXM
figure('Name', 'Fit Ct: 2CXM');
subplot(2,2,1);
plot(F_p, model_params_C_2CXM(:,1), 'rx');
title('F_p');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,2);
plot(PS, model_params_C_2CXM(:,2), 'rx');
title('PS');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,3);
plot(v_e, model_params_C_2CXM(:,3), 'rx');
title('v_e');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,4);
plot(v_p, model_params_C_2CXM(:,4), 'rx');
title('v_p');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

%Look at Ve and Vp estimates for ETM
figure('Name', 'Fit Ct: ETM');
subplot(1,2,1);
plot(v_e, model_params_C_ETM(:,2), 'rx');
title('v_e');
xlabel('Actual');
ylabel('Predicted (ETM)');
axis equal

subplot(1,2,2);
plot(v_p, model_params_C_ETM(:,3), 'rx');
title('v_p');
xlabel('Actual');
ylabel('Predicted (ETM)');
axis equal
%%
% 3) Fit signals using actual T1
[model_params_S, model_fit_S0, ~,~,CmS_t] = run_madym_lite(...
    '2CXM', S_tn, 'dyn_times', dyn_t,...
    'input_Ct', 0,...
    'T1', T1,...
    'TR', TR,...
    'FA', dyn_FA,...
    'r1Const', r1Const,...
    'injectionImage', injection_img);

%Convert the modelled concentrations back to signal space
Sm_t = concentration_to_signal(...
    CmS_t, dyn_FA, TR, T1, S0, r1Const, false);

display_dce_model_fits(S_tn, Sm_t, ...
    'plot_rows', 4, 'plot_cols', 6);
%
% Look at parameter estimates for 2CXM
figure('Name', 'Fit St: 2CXM, known T1');
subplot(2,2,1);
plot(F_p, model_params_S(:,1), 'rx');
title('F_p');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,2);
plot(PS, model_params_S(:,2), 'rx');
title('PS');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,3);
plot(v_e, model_params_S(:,3), 'rx');
title('v_e');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,4);
plot(v_p, model_params_S(:,4), 'rx');
title('v_p');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal
%%
% 4) Fit signals using estimated T1
[model_params_S, model_fit_S1, ~,~,CmS_t] = run_madym_lite(...
    '2CXM', S_tn, 'dyn_times', dyn_t,...
    'input_Ct', 0,...
    'T1', T1_fit,...
    'TR', TR,...
    'FA', dyn_FA,...
    'r1Const', r1Const,...
    'injectionImage', injection_img);

%Convert the modelled concentrations back to signal space
Sm_t = concentration_to_signal(...
    CmS_t, dyn_FA, TR, T1_fit, S0_fit, r1Const, false);

display_dce_model_fits(S_tn, Sm_t, ...
    'plot_rows', 4, 'plot_cols', 6);
%
% Look at parameter estimates for 2CXM
figure('Name', 'Fit St: 2CXM, estimated T1');
subplot(2,2,1);
plot(F_p, model_params_S(:,1), 'rx');
title('F_p');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,2);
plot(PS, model_params_S(:,2), 'rx');
title('PS');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,3);
plot(v_e, model_params_S(:,3), 'rx');
title('v_e');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,4);
plot(v_p, model_params_S(:,4), 'rx');
title('v_p');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal
%%
% 5) Fit signals using estimated T1 and S0
[model_params_S, model_fit_S2, ~,~,CmS_t] = run_madym_lite(...
    '2CXM', S_tn, 'dyn_times', dyn_t,...
    'input_Ct', 0,...
    'T1', T1_fit,...
    'S0', S0_fit,...
    'useRatio', 0,...
    'TR', TR,...
    'FA', dyn_FA,...
    'r1Const', r1Const,...
    'injectionImage', injection_img,...
    'dummy_run', 0);

%Convert the modelled concentrations back to signal space
Sm_t = concentration_to_signal(...
    CmS_t, dyn_FA, TR, T1_fit, S0_fit, r1Const, false);

display_dce_model_fits(S_tn, Sm_t, ...
    'plot_rows', 4, 'plot_cols', 6);
%
% Look at parameter estimates for 2CXM
figure('Name', 'Fit St: 2CXM, estimated T1');
subplot(2,2,1);
plot(F_p, model_params_S(:,1), 'rx');
title('F_p');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,2);
plot(PS, model_params_S(:,2), 'rx');
title('PS');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,3);
plot(v_e, model_params_S(:,3), 'rx');
title('v_e');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

subplot(2,2,4);
plot(v_p, model_params_S(:,4), 'rx');
title('v_p');
xlabel('Actual');
ylabel('Predicted (2CXM)');
axis equal

%%
fprintf('Mean model fit residual, concentration, 2CXM: %5.4f\n',...
    mean(model_fit_C_2CXM));
fprintf('Mean model fit residual, concentration, ETM: %5.4f\n',...
    mean(model_fit_C_ETM));
fprintf('Mean model fit residual, signal, 2CXM, known T1: %5.4f\n',...
    mean(model_fit_S0));
fprintf('Mean model fit residual, signal, 2CXM, estimated T1, ratio method: %5.4f\n',...
    mean(model_fit_S1));
fprintf('Mean model fit residual, signal, 2CXM, estimated T1, S0: %5.4f\n',...
    mean(model_fit_S2));