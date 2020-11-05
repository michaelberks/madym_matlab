figure;
subplot(2, 3, 1); 
plot(model_params(:,1,3), model_params(:,2,3), 'r.');
xlabel('F+'); ylabel('F-'); 

subplot(2, 3, 2); 
plot(model_params(:,1,3), model_params(:,3,3), 'r.');
xlabel('F+'); ylabel('K+'); 

subplot(2, 3, 3); 
plot(model_params(:,1,3), model_params(:,4,3), 'r.');
xlabel('F+'); ylabel('K-'); 

subplot(2, 3, 4); 
plot(model_params(:,2,3), model_params(:,3,3), 'r.');
xlabel('F-'); ylabel('K+'); 

subplot(2, 3, 5); 
plot(model_params(:,2,3), model_params(:,4,3), 'r.');
xlabel('F-'); ylabel('K-'); 

subplot(2, 3, 6); 
plot(model_params(:,3,3), model_params(:,4,3), 'r.');
xlabel('K+'); ylabel('K-'); 
%%
milano_monte_carlo_fixed(...
    'experiment_root', experiment_root,...
    'results_dir', 'fixed',...
    'analysis_dir', 'analysis_irf2',...
    'gt_idx', [1 3],...
    'irf2_2cxm', 1,...
    'n_pts', 1e3,...
    'n_t', 0,...
    'temporal_resolution', 3.8,...
    'end_time', 7,...
    'injection_image', 0,...
    'injection_time', 60,...
    'do_fitting', 0,...
    'do_analysis', 1,...
    'dummy_run', 0);
%%
v_ecs = v_ecs_active;
k_i = k_i_active;
k_ef = linspace(0, 2*k_i_active, 20)'; 

[Fpos, Fneg, Kpos, Kneg] = ...
    active_params_phys_to_model(...
    F_p_active, v_ecs, k_i, k_ef(:));

[F_p, PS, v_e, v_p] = ...
    two_cx_params_model_to_phys(Fpos, Fneg, Kpos, Kneg);
%%
Eneg = Fneg ./ F_p;
Tb = 1 ./ (Kpos - Eneg.*(Kpos-Kneg));
Tei = Tb .* Kpos .* Kneg;
Tp = 1 ./ (Kpos + Kneg - Tei);

v = F_p.*(Tb - Tp) ./ ...
    (Tp.*Tb.*Kpos.*Kneg);
%%
v_e2 = F_p.*(Tb.*(Kpos + Kneg - Tb.*Kpos.*Kneg) - 1) ./ ...
    (Tb.*Kpos.*Kneg);
%%
Tbi = 1 ./ Tb;
v_e3 = F_p.*(Kpos + Kneg - Tb.*Kpos.*Kneg - Tbi) ./ ...
    (Kpos.*Kneg);
%%
v_e4 = F_p .* (Kneg - Kpos.*Kneg ./ (Kpos - Eneg.*(Kpos-Kneg)) + Eneg.*(Kpos-Kneg)) ./...
    (Kpos.*Kneg);
%%
v_e5 = F_p .* (Kneg - Kpos.*Kneg.*(Fpos + Fneg)./(Kpos.*Fpos + Kneg.*Fneg) + ...
    Fneg.*(Kpos - Kneg)./(Fpos + Fneg)) ./...
    (Kpos.*Kneg);
%%
v_e6 = (Fpos+Fneg)./Kpos - (Fpos+Fneg).^2 ./ (Kpos.*Fpos + Kneg.*Fneg) +...
    Fneg.*(Kpos-Kneg)./(Kpos.*Kneg);
%%
Epos = Fpos ./ F_p;
Ti = (1 - v_ecs)./k_ef;
Te = v_ecs ./ (F_p + k_i);
v_e7 = F_p.*(Ti - 1./(Epos./Ti + (1-Epos)./Te) + (1-Epos).*(Te-Ti));
%%
v_e8 = F_p.*(Ti - 1./(Epos./Ti + (1-Epos)./Te) + (1-Epos).*(Te-Ti));
%%
v_e9 = F_p.*(Ti - Ti.*Te./(Te.*Epos + Ti.*(1-Epos)) + (1-Epos).*(Te-Ti));
%%
v_e10 = F_p.*(Ti - Te./(1-E_i) + (1-Epos).*(Te-Ti));
%%
v_e11 = F_p.*(Te - Te./(1-E_i) + E_i.*Ti);
%%
v_e12 = F_p.*(Te.*E_i./(E_i-1) + E_i.*Ti);
%%
v_e13 = F_p.*E_i.*(Te./(E_i-1) + Ti);
%%
v_e14 = F_p.*k_i.*((1-v_ecs)./k_ef - v_ecs./F_p)./(F_p+k_i);
%%
v_p2 = F_p ./ (Kneg + Epos.*(Kpos-Kneg));
%%
v_p3 = F_p ./ (1./Te + Epos.*(1./Ti - 1./Te));
%%
v_p4 = F_p ./ (1./Te + E_i.*(1./Ti - 1./Te)./(1 - Te./Ti));
%%
v_p5 = Ti.*F_p ./ (Ti./Te + E_i.*(Te-Ti)./((Ti-Te).*Te.*Ti));
%%
v_p6 = F_p.*Te./(1-E_i);
%%
PS2 = F_p.*(Tb.*Kpos + Tb.*Kneg - Tb.*Tb.*Kpos.*Kneg - 1);
%%
PS3 = Tb.*F_p.*(Kpos - Tb.*Kpos.*Kneg - Epos.*(Kpos-Kneg));
%%
PS4 = F_p.*(Kpos - Kpos.*Kneg./(Epos.*(Kpos-Kneg)+Kneg) - Epos.*(Kpos-Kneg))./(Epos.*(Kpos-Kneg)+Kneg);
%%
PS5 = F_p.*Te.*(1./Ti - 1./(Ti.*(1-Ei)) + Ei./Te)./(1- Ei);
%%
PS6 = (k_i+F_p).*(Te./Ti - Te.*(k_i+F_p)./(Ti.*F_p) + k_i./(k_i + F_p));
%%
PS7 = (k_i+F_p).*((v_ecs.*k_ef)./((k_i+F_p).*(1-v_ecs)) - v_ecs.*k_ef./((1-v_ecs).*F_p) + k_i./(k_i + F_p));
%%
PS8 = k_i .* (1 - v_ecs.*k_ef./(1-v_ecs)./F_p);
