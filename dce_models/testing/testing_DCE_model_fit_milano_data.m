visit_dir = 'C:\isbe\qbi\data\milano_primovist\PRIMDCE_8\visit1\';
etm_dir = [visit_dir 'mdm_analysis_reg_VP\'];
mam_dir = [visit_dir 'mdm_analysis_reg_MATERNE\'];
gsm_dir = [visit_dir 'mdm_analysis_reg_GADOXETATE\'];
%%
etm.ktrans = load_img_volume([etm_dir 'Ktrans.hdr']);
etm.ve = load_img_volume([etm_dir 'Ve.hdr']);
etm.vp = load_img_volume([etm_dir 'Vp.hdr']);
etm.offset = load_img_volume([etm_dir 'offset.hdr']);
etm.err = load_img_volume([etm_dir 'ERR.hdr']);
etm.error_img = load_img_volume([etm_dir 'error_img.hdr']);
%%
mam.k1a = load_img_volume([mam_dir 'k1a.hdr']);
mam.k1v = load_img_volume([mam_dir 'k1hpv.hdr']);
mam.k2 = load_img_volume([mam_dir 'k2.hdr']);
mam.aoffset = load_img_volume([mam_dir 'aoffset.hdr']);
mam.voffset = load_img_volume([mam_dir 'voffset.hdr']);
mam.err = load_img_volume([mam_dir 'ERR.hdr']);
mam.error_img = load_img_volume([mam_dir 'error_img.hdr']);
%%
gsm.Fp = load_img_volume([gsm_dir 'Fp.hdr']);
gsm.fa = load_img_volume([gsm_dir 'fa.hdr']);
gsm.vecs = load_img_volume([gsm_dir 'vecs.hdr']);
gsm.kin = load_img_volume([gsm_dir 'kin.hdr']);
gsm.kout = load_img_volume([gsm_dir 'kout.hdr']);
gsm.aoffset = load_img_volume([gsm_dir 'aoffset.hdr']);
gsm.voffset = load_img_volume([gsm_dir 'voffset.hdr']);
gsm.err = load_img_volume([gsm_dir 'ERR.hdr']);
%%
roi = load_img_volume([visit_dir 'liver_top_mask.hdr']) > 0;
roi_liver = roi;
%%
aif = load_aif([visit_dir 'mdm_analysis_T1_ants\slice_19_Auto_AIF.txt']);
T1 = load_img_volume([visit_dir 'mdm_analysis_T1_ants\T1.hdr']);

n_vols = size(aif, 1);
[dyn_times, dyn_FA, dyn_TR] = get_dyn_times([visit_dir 'dynamic\dyn_'], n_vols);
aif = aif(:,2);
[pif] = compute_PIF(aif, [], dyn_times);
aif_injection_image = 8;
relax_coeff = 14.0e-3;
Hct = 0.42;
%%
dyn_signals = get_dyn_vals([visit_dir 'dynamic_reg_ants\dyn_'], n_vols, roi);
dyn_concentrations = ...
    signal_to_concentration(dyn_signals, dyn_FA(1), dyn_TR(1), ...
    T1(roi),...
    relax_coeff, aif_injection_image);
%%
[Cl_etm] =...
    extended_kety_model(dyn_times, aif / (1 - Hct),...
    etm.ktrans(roi),...
    etm.vp(roi),...
    etm.ve(roi),...
    etm.offset(roi));
etm.err_mat = sum((dyn_concentrations - Cl_etm).^2);

[Cl_gsm] = gadoxetate_model(...
    gsm.Fp(roi)',...
    gsm.vecs(roi)',...
    gsm.kin(roi)',...
    gsm.kout(roi)',...
    gsm.fa(roi)',...
    gsm.aoffset(roi)',...
    aif, pif, dyn_times, Hct);
gsm.err_mat = sum((dyn_concentrations - Cl_gsm).^2);
%%
[Cl_mam] = ...
    materne_model(...
    mam.k1a(roi),...
    mam.k1v(roi),...
    mam.k2(roi),...
    mam.aoffset(roi),...
    aif, pif, dyn_times, Hct, mam.voffset(roi));
mam.err_mat = sum((dyn_concentrations - Cl_mam).^2);
%%
%As a sanity check compare the model fit errors we generate in matlab vs 
%the errors returned by madym 
figure;
subplot(2,3,1);
plot(etm.err(roi), etm.err_mat, 'r.'); axis equal;
title('Extended-Tofts model');
xlabel('Madym SSE');
ylabel('Matlab SSE');
axis([0 1.2 0 1.2]);

subplot(2,3,2);
plot(gsm.err(roi), gsm.err_mat, 'r.'); axis equal;
title('Gadoxetate specific model');
xlabel('Madym SSE');
ylabel('Matlab SSE');
axis([0 1.2 0 1.2]);

subplot(2,3,3);
plot(mam.err(roi), mam.err_mat, 'r.'); axis equal;
title('Materne model');
xlabel('Madym SSE');
ylabel('Matlab SSE');
axis([0 1.2 0 1.2]);

%Compare the model fit errors between each other
subplot(2,3,4);
plot(etm.err_mat, gsm.err_mat, 'r.'); axis equal; hold on;
plot([0 1.2], [0 1.2], 'k--');
title('ETM vs GSM');
xlabel('ETM SSE');
ylabel('GSM SSE');
axis([0 1.2 0 1.2]);

subplot(2,3,5);
plot(etm.err_mat, mam.err_mat, 'r.'); axis equal; hold on;
plot([0 1.2], [0 1.2], 'k--');
title('ETM vs Materne');
xlabel('ETM SSE');
ylabel('Materne SSE');
axis([0 1.2 0 1.2]);

subplot(2,3,6);
plot(gsm.err_mat, mam.err_mat, 'r.'); axis equal; hold on;
plot([0 1.2], [0 1.2], 'k--');
title('GSM vs Materne');
xlabel('GSM SSE');
ylabel('Materne SSE');
axis([0 1.2 0 1.2]);
%%
T1_roi = T1(roi);
%%
[k1a, k1v, k2, aoffset, voffset, model_signals] = fit_materne_model(...
    dyn_signals(:,1), dyn_times, aif/(1-Hct), pif/(1-Hct), T1_roi(1),...
    'FA', dyn_FA(1), 'TR', dyn_TR(1), 'num_itr', 100, ...
    'num_baseline', aif_injection_image, 'relax_coeff', relax_coeff,...
    'fit_to', 'concentration');
[Cl_mam1] = ...
        materne_model(k1a, k1v, k2, aoffset, aif/(1-Hct), pif/(1-Hct), dyn_times, 0, voffset);
%%    
figure; 
a1 = subplot(1,2,1);
plot(dyn_times, dyn_concentrations(:,1)); hold all;
plot(dyn_times, Cl_mam1(:,1), 'k', 'linewidth', 3);
title(sprintf('SSE = %f', sum((dyn_concentrations(:,1)-Cl_mam1).^2)));

a2 = subplot(1,2,2);
plot(dyn_times, dyn_signals(:,1)); hold all;
plot(dyn_times, model_signals(:,1), 'k', 'linewidth', 3);
title(sprintf('SSE = %f', sum((dyn_signals(:,1)-model_signals).^2)));

sse = zeros(21, 1); idx = 1;
for k1v = linspace(0, 1, 21)
    [Cl_mam1] = ...
        materne_model(k1a, k1v, k2, aoffset, aif/(1-Hct), pif/(1-Hct), dyn_times, 0, voffset);


    plot(a1, dyn_times, Cl_mam1(:,1), '--');
    
    sse(idx, 1) = sum((dyn_concentrations(:,1)-Cl_mam1).^2);
    idx = idx + 1;
end

%%
[k1a_fit, k1v_fit, k2_fit, aoffset_fit, voffset_fit] = fit_materne_model(...
    dyn_signals(:,1:30), dyn_times, aif/(1-Hct), pif/(1-Hct), T1_roi(1:30),...
    'FA', dyn_FA(1), 'TR', dyn_TR(1), 'num_itr', 100, ...
    'num_baseline', aif_injection_image, 'relax_coeff', relax_coeff,...
    'fit_to', 'concentration');

[Cl_mam_fit] = ...
    materne_model(k1a_fit, k1v_fit, k2_fit,...
    aoffset_fit, aif/(1-Hct), pif/(1-Hct), dyn_times, 0, voffset_fit);

display_dce_model_fits(...
    dyn_concentrations(:,1:30), Cl_mam_fit);
%%
figure;
subplot(1,2,1);
imgray(mam.k1a(:,:,18)); colorbar;
subplot(1,2,2);
imgray(mam.k1a(1:64,1:80,18)); colorbar; caxis([0 0.5]);
title('k_{1a}');
colormap jet;
%
figure;
subplot(1,2,1);
imgray(mam.k1v(:,:,18)); colorbar;
subplot(1,2,2);
imgray(mam.k1v(1:64,1:80,18)); colorbar; caxis([0 0.5]);
title('k_{1hpv}');
colormap jet;
%
figure;
subplot(1,2,1);
imgray(mam.k2(:,:,18)); colorbar;
subplot(1,2,2);
imgray(mam.k2(1:64,1:80,18)); colorbar; caxis([0 1.0]);
title('k_2');
colormap jet;
%
figure;
subplot(1,2,1);
imgray(mam.aoffset(:,:,18)); colorbar;
subplot(1,2,2);
imgray(mam.aoffset(1:64,1:80,18)); colorbar; caxis([0 0.5]);
title('t_a');
colormap jet;
%
figure;
subplot(1,2,1);
imgray(mam.voffset(:,:,18)); colorbar;
subplot(1,2,2);
imgray(mam.voffset(1:64,1:80,18)); colorbar; caxis([0 0.5]);
title('t_v');
colormap jet;
%
%%
for i_sub = 1:10
    visit_dir = ['C:\isbe\qbi\data\milano_primovist\PRIMDCE_' num2str(i_sub) '\visit1\'];
    mam_dir = [visit_dir 'mdm_analysis_reg_MATERNE\'];
    voffset = load_img_volume([mam_dir 'voffset.hdr']);
    
    figure;
    subplot(1,2,1);
    imgray(voffset(:,:,18)); colorbar;
    subplot(1,2,2);
    imgray(voffset(1:64,1:80,18)); colorbar; caxis([0 0.5]);
    title('t_v');
    colormap jet;
end
%%
[roi_analysis] = milano_analysis(8, ... 
    'reg', 1,...
    'do_save', 0,...
    'make_tumour_roi', [],...
    'make_liver_roi', [],...
    'do_region_analysis', 0,...
    'do_point_selection', 0,...
    'do_tse_figs', 0,...
    'do_aic_bin_figs', 1,...
    'do_etm_maps', 0,...
    'do_gsm_maps', 0,...
    'do_mam_maps', 1,...
    'do_roi_overlay_figs', 0, ...
    'do_series_figs', 0,...
    'do_mask_figs', 0,...
    'do_mask_comp_figs', 0,...
    'do_aic_figs', 0,...
    'do_etm_overlay', 0,...
    'do_gsm_overlay', 0,...
    'Hct', 0.42);
%%
mam2_dir = [visit_dir 'fri_test_mat2_reg_MATERNE\'];
mam2.k1a = load_img_volume([mam2_dir 'k1a.hdr']);
mam2.k1v = load_img_volume([mam2_dir 'k1hpv.hdr']);
mam2.k2 = load_img_volume([mam2_dir 'k2.hdr']);
mam2.aoffset = load_img_volume([mam2_dir 'aoffset.hdr']);
mam2.voffset = load_img_volume([mam2_dir 'voffset.hdr']);
mam2.err = load_img_volume([mam2_dir 'ERR.hdr']);
mam2.error_img = load_img_volume([mam2_dir 'error_img.hdr']);
%%
[Cl_mam_k1a] = ...
    materne_model(...
    mam.k1a(roi),...
    mam.k1v(roi),...
    mam.k2(roi),...
    mam.aoffset(roi),...
    aif, pif, dyn_times, Hct, mam.voffset(roi)+1e-6);
err_mat_k1a = sum((dyn_concentrations - Cl_mam_k1a).^2);
mean(abs(err_mat_k1a - mam.err_mat))