%Check DIIRF
% The 2CXM and Gadoxetate models should both be reproduced by DIIRF

%Set up times and IFs
n_t = 100;
t = linspace(0,6,n_t);
Ca_t = population_aif(t, 8);
Cv_t = compute_PIF(Ca_t, [], t);

%%
%Shared params
F_p = 1.0;

%2CXM params
v_e = 0.2;
v_p = 0.1;
PS = 0.3;
[K_pos_2cxm, K_neg_2cxm, F_pos_2cxm, F_neg_2cxm] = ...
    two_cx_params_phys_to_model(F_p, PS, v_e, v_p, false);

%Gadoxetate params
v_ecs = 0.2;
k_i = 0.3;
k_ef = 0.01;
f_a = 0.3;
[K_pos_gsm, K_neg_gsm, F_pos_gsm, F_neg_gsm] = ...
    active_params_phys_to_model(F_p, v_ecs, k_i, k_ef);

C_gsm = gadoxetate_model(...
    F_p, v_ecs, k_i, k_ef, f_a, 0, Ca_t, Cv_t, t, 0);
C_2cxm = two_cx_model(...
    F_p, PS, v_e, v_p,  0, Ca_t, t);

%Convert active uptake physiological parameters to their functional form
%then apply the DI-IRF model, this should be identical to the gadoxetate
%model
C_gsm_diirf = diirf_model(...
    F_pos_gsm, F_neg_gsm, K_pos_gsm, K_neg_gsm, f_a, 0, Ca_t, Cv_t, t, 0);

%Convert exchange physiological parameters to their functional form
%then apply the DI-IRF model, this should be identical to the
%two-compartment exchange model. Note we need to set f_a = 1.0, to match
%the single input of the 2CXM
C_2cxm_diirf = diirf_model(...
    F_pos_2cxm, F_neg_2cxm, K_pos_2cxm, K_neg_2cxm, 1.0, 0, Ca_t, Cv_t, t, 0);
%%
figure;
subplot(1,2,1);
title('Gadoxetate vs DI-IRF, f_a = 1'); hold all;
plot(t, C_gsm, 'b', 'linewidth', 3);
plot(t, C_gsm_diirf, 'r--', 'linewidth', 3.0);
legend({'GSM', 'DI-IRF'});

subplot(1,2,2);
title('2CXM vs DI-IRF , f_a = 1'); hold all;
plot(t, C_2cxm, 'b', 'linewidth', 3);
plot(t, C_2cxm_diirf, 'r--', 'linewidth', 3.0);
legend({'2CXM', 'DI-IRF'});