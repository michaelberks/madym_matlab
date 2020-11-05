%For single input models we start with the 2CXM model as the model general,
%setting F_p = Inf should produce the ETM, and setting vp = 0 produce the
%TM. Similarly setting vp = 0 in the ETM should produce the TM.
%
% In the dual input models, starting with the most general, Leo's full
% gadoxetate model, setting k_ef = 0 produces the model described in the
% Sourbron 2012 radiology paper, and further setting k_i = 0 reduces to the
% single compartment Materne model. We sanity check these results below
%%
%% First set up some time-points and input functions
n_t = 100;
t = linspace(0,6,n_t);
Ca_t = population_aif(t, 8);
Cv_t = compute_PIF(Ca_t, [], t);
%% Check the single input models: 2CXM vs ETM

%% 1) When vp = 0, all other params positive finite
% 2CXM reduces to TM, ve = ve, Ktrans = (F_p.*PS) ./ (F_p + PS)
F_p = 2.0;
PS = 0.20;
v_e = 0.3;
v_p = 0.0;
C_2cxm = two_cx_model(...
    F_p, PS, v_e, v_p,  0, Ca_t, t);
%
Ktrans = (F_p.*PS) ./ (F_p + PS);
C_etm = extended_tofts_model(Ktrans, v_e, v_p, 0, Ca_t, t);
   
figure;
title('ETM vs 2CXM v_p = 0, 0 < F_p,v_e,PS < \infty'); hold all;
plot(t, C_2cxm, 'linewidth', 3.0);
plot(t, C_etm, '--', 'linewidth', 2.0);
legend({'2CXM', 'TM'});
%% 2) When ve = 0, all other params positive finite
% 2CXM reduces to TM, with ve = vp, Ktrans = F_p, vp = 0
F_p = 0.823;
PS = 0.20;
v_e = 0.0;
v_p = 0.3;
C_2cxm = two_cx_model(...
    F_p, PS, v_e, v_p,  0, Ca_t, t);
%
Ktrans = F_p;
v_e_etm = v_p;
v_p_etm = 0;

C_etm = extended_tofts_model(Ktrans, v_e_etm, v_p_etm, 0, Ca_t, t);

figure;
title('ETM vs 2CXM v_e = 0, 0 < F_p,v_e,PS < \infty'); hold all;
plot(t, C_2cxm, 'linewidth', 3.0);
plot(t, C_etm, '--', 'linewidth', 2.0);
legend({'2CXM', 'TM'});
%% 3) When PS = 0, all other params positive finite
% 2CXM reduces to TM, with ve = vp, Ktrans = F_p, vp = 0
F_p = 0.4/(1-0.42);
PS = 0.0;
v_e = 0.3;
v_p = 0.032;
C_2cxm = two_cx_model(...
    F_p, PS, v_e, v_p,  0, Ca_t, t);
%
Ktrans = F_p;
v_e_etm = v_p;
v_p_etm = 0;

C_etm = extended_tofts_model(Ktrans, v_e_etm, v_p_etm, 0, Ca_t, t);
  
figure;
title('ETM vs 2CXM PS = 0, 0 < F_p,v_e,v_p < \infty'); hold all;
plot(t, C_2cxm, 'linewidth', 3.0);
plot(t, C_etm, '--', 'linewidth', 2.0);
legend({'2CXM', 'TM'});
%% 4) When Fp -> Inf, all other params positive finite
% 2CXM reduces to ETM, with ve = ve, vp = vp, Ktrans = PS
F_p = [1 2 3 4 5 6];
PS = 0.20;
v_e = 0.4;
v_p = 0.2;
C_2cxm = two_cx_model(...
    F_p, PS, v_e, v_p,  0, Ca_t, t);

Ktrans = PS;
v_e_etm = v_e;
v_p_etm = v_p;

C_etm = extended_tofts_model(Ktrans, v_e_etm, v_p_etm, 0, Ca_t, t);
   
figure;
title('ETM vs 2CXM F_p \rightarrow \infty, 0 < v_p,v_e,PS < \infty'); hold all;
plot(t, C_etm, 'k-', 'linewidth', 3.0);
plot(t, C_2cxm, '--', 'linewidth', 1.0);
labels = cell(length(F_p) + 1, 1);
labels{1} = 'ETM';
for i_fp = 1:length(F_p)
    labels{i_fp+1} = sprintf('2CXM, F_p = %3.2f', F_p(i_fp));
end
legend(labels);

%% 5) When PS -> Inf, all other params positive finite
% 2CXM reduces to ETM, with ve = ve + vp, vp = vp, Ktrans = FP
F_p = 1.0;
PS = 1e9;
v_e = 0.3;
v_p = 0.2;
C_2cxm = two_cx_model(...
    F_p, PS, v_e, v_p,  0, Ca_t, t);
%
Ktrans = F_p;
v_e_etm = v_e + v_p;
v_p_etm = 0;

C_etm = extended_tofts_model(Ktrans, v_e_etm, v_p_etm, 0, Ca_t, t);
  
figure;
title('ETM vs 2CXM PS \rightarrow \infty, 0 < v_p,v_e,FP < \infty'); hold all;
plot(t, C_2cxm, 'linewidth', 3.0);
plot(t, C_etm, '--', 'linewidth', 2.0);
legend({'2CXM', 'ETM'});
%% Check the dual-input models: Gadoxetate vs Materne
F_p = 1.0;
v_e = 0.2;
k_i = [0.0 0.001 0.1];
k_ef = 0;
f_a = 0.5;

%These parameters allow us to derive the inputs for the Materne model
k1a = F_p*f_a;
k1v = F_p*(1-f_a);
k2 = F_p / v_e;

C_mm = materne_model(k1a, k1v, k2, 0, Ca_t, Cv_t, t, 0);
C_gsm = gadoxetate_model(F_p, v_e, k_i, k_ef, f_a, 0, Ca_t, Cv_t, t, 0);

figure;
title('Gadoxetate vs Materne, K_i = K_{ef} = 0, K_2 = F_p / v_e'); hold all;
plot(t, C_gsm, 'linewidth', 3);
plot(t, C_mm, 'k--', 'linewidth', 2.0);
legend({'GSM K_i = 0.0', 'GSM K_i = 0.001', 'GSM K_i = 0.1', 'Materne'});
%%
F_p = [0.1 0.5 1.0];
v_e = 1.0;
k_i = 0.01;
k_ef = 0.02;
f_a = 0.7;

%These parameters allow us to derive the inputs for the Materne model
k1a = F_p*f_a;
k1v = F_p*(1-f_a);
k2 = F_p / v_e;

C_mm = materne_model(k1a, k1v, k2, 0, Ca_t, Cv_t, t, 0);
C_gsm = gadoxetate_model(F_p, v_e, k_i, k_ef, f_a, 0, Ca_t, Cv_t, t, 0);

figure;
title('Gadoxetate vs Materne, K_i = K_{ef} = 0.01, K_2 = F_p / v_e'); hold all;
plot(t, C_gsm, 'linewidth', 3);
plot(t, C_mm, 'k--', 'linewidth', 2.0);
legend({'GSM F_p = 0.1', 'GSM K_p = 0.5', 'GSM F_p = 1.0', 'Materne'});

%% Check single vs dual: Materne vs TM
Ktrans = 0.25;
v_e = 0.5;

%Materne params are derived as below, note we need make sure we turn off
%Hct correction in the materne model, as we already have this in our AIF;
k1a = Ktrans;
k1v = 0;
k2 = Ktrans / v_e;

C_tm = extended_tofts_model(Ktrans, v_e, 0, 0, Ca_t, t);
C_mm = materne_model(k1a, k1v, k2, 0, Ca_t, Cv_t, t, 0);

figure; hold all;
title('TM vs Materne f_a = 1, F_p = K_{trans}, K_2 = K_{trans} / v_e');
plot(t, C_tm, 'linewidth', 3);
plot(t, C_mm, 'r--', 'linewidth', 2.0);
legend({'TM', 'Materne'});

%% Check single vs dual: Gadoxetate with f_a = 1, k_ef = 0 vs ETM ????
F_p = 2.0;
v_e = 0.2;
k_i = 0.0;
k_ef = 0;
f_a = 1.0;

%These parameters allow us to derive the inputs for the ETM
E_i = k_i ./ (F_p + k_i);
T_e = v_e ./ (F_p + k_i);

v_p = F_p.*E_i; %if k_i = 0, v_p = 0
v_etm = F_p.*(1 - E_i) .* T_e; %if k_i = 0, v_etm = v_e .* F_p
Ktrans = F_p.*(1 - E_i); %if k_i = 0, ktrans = 1.0

C_etm = extended_tofts_model(Ktrans, v_e, v_p, 0, Ca_t, t);
C_gsm = gadoxetate_model(F_p, v_e, k_i, k_ef, f_a, 0, Ca_t, Cv_t, t, 0);

figure;
title('Gadoxetate vs ETM, f_a = 1'); hold all;
plot(t, C_gsm, 'linewidth', 3);
plot(t, C_etm, '--', 'linewidth', 3.0);
legend({'GSM', 'ETM'});

%% Check single vs dual: Gadoxetate with f_a = 1 vs ETM
F_p = 2.0;
v_e = 0.2;
k_i = 0.3;
k_ef = [0.4 0.31 0.30001 0.3 0.2999999 0.2];
f_a = 1.0;

E_i = k_i ./ (F_p + k_i);
T_e = v_e ./ (F_p + k_i);
T_i = (1 - v_e) ./ k_ef;

E_neg = E_i ./ (1 - T_e ./ T_i);
F_neg = F_p .* E_neg;
F_pos = F_p .* (1 - E_neg);
K_neg = 1 ./ T_i;
K_pos = 1 ./ T_e;

[F_p2, PS2, v_e2, v_p2] = ...
    two_cx_params_model_to_phys(K_pos, K_neg, F_pos, F_neg, false);

C_gsm = gadoxetate_model(...
    F_p, v_e, k_i, k_ef, f_a, 0, Ca_t, Cv_t, t, 0);
C_2cxm = two_cx_model(...
    F_p2, PS2, v_e2, v_p2,  0, Ca_t, t);

figure;
title('Gadoxetate vs 2CXM, f_a = 1'); hold all;
plot(t, C_gsm, 'b', 'linewidth', 3);
plot(t, C_2cxm, 'r--', 'linewidth', 3.0);
legend({'GSM', '2CXM'});
%%
F_p = 2.0;
v_e = 0.2;
k_i = 0.3;
k_ef = 10.^(-1*(4:8));
f_a = 1.0;

E_i = k_i ./ (F_p + k_i);
T_e = v_e ./ (F_p + k_i);
T_i = (1 - v_e) ./ k_ef;

E_neg = E_i ./ (1 - T_e ./ T_i);
F_neg = F_p .* E_neg;
F_pos = F_p .* (1 - E_neg);
K_neg = 1 ./ T_i;
K_pos = 1 ./ T_e;

figure; plot(k_ef, F_neg ./ K_neg); hold on;
plot(k_ef, F_neg ./ K_neg, 'rx');



    