% Showing that the active model derivation is only valid for beta_pos < beta_neg,
% otherwise K_i < 0. In this scenario we can swap the 'positive' and
% 'negative' components of the bi-exponential, generating the same
% concentration curve, but with a derivation where K_i > 0;
%%
alpha_pos = 4*rand(1e4,1);
alpha_neg = 4*rand(1e4,1);
beta_pos = 4*rand(1e4,1);
beta_neg = 4*rand(1e4,1);
%%
%Apply the active model derivations
F_p = alpha_pos + alpha_neg;
E_pos = alpha_pos ./ F_p;
T_e = 1 ./ beta_neg;
T_i = 1 ./ beta_pos;
E_i = E_pos.*(1 - T_e./T_i);

K_i = E_i.*F_p ./ (1 - E_i);
v_ecs = F_p.*T_e ./ (1 - E_i);
%%
% Let a*_+ = a_-, a*_- = a_+, b*_+ = b*_-, b*_- = b_+ be the 'swapped' 
% derivation
%
%It's (moderately) easy to show F_p2 = F_p and v_ecs2 = v_ecs from the
%active uptake equations, confirmed numerically below
alpha_pos2 = alpha_neg;
alpha_neg2 = alpha_pos;
beta_pos2 = beta_neg;
beta_neg2 = beta_pos;
%
F_p2 = alpha_pos2 + alpha_neg2;
E_pos2 = alpha_pos2 ./ F_p2;
T_e2 = 1 ./ beta_neg2;
T_i2 = 1 ./ beta_pos2;
E_i2 = E_pos2.*(1 - T_e2./T_i2);

K_i2 = E_i2.*F_p2 ./ (1 - E_i2);
v_ecs2 = F_p2.*T_e2 ./ (1 - E_i2);

figure;
subplot(1,2,1);
plot(F_p, F_p2, 'r.'); axis equal;
xlabel('F_p');
ylabel('F_p^*');
title('Derived F_p when we set \alpha,\beta_{+/-} = \alpha,\beta_{-/+}')
axis([0 8 0 8]);

subplot(1,2,2);
plot(v_ecs, v_ecs2, 'r.'); axis equal;
xlabel('v_{ecs}');
ylabel('v_{ecs}^*');
title('Derived v_{ecs} when we set \alpha,\beta_{+/-} = \alpha,\beta_{-/+}')
axis([0 100 0 100]);

%%
%We will show K_i > 0 iff K_i2 < 0, as can be seen numerically
figure;
plot(K_i, K_i2, 'r.'); axis equal;
xlabel('K_i');
ylabel('K_i^*');
title('Derived K_i when we set \alpha,\beta_{+/-} = \alpha,\beta_{-/+}')
axis([-10 60 -10 60]);
%%
%Firstly, from the active derivation equations we show 
% K_i > 0 iff b_+ < b_-, which can also be seen numerically
figure;
plot(K_i, beta_neg - beta_pos, 'r.');
xlabel('K_i');
ylabel('\beta_- - \beta_+');
title('Derived K_i when we set \alpha,\beta_{+/-} = \alpha,\beta_{-/+}')
axis([-10 60 -4 4]);
%Hence it suffices to show K_i2 < 0 iff b_+ < b_-

%Next we derive K_i2 = A*K_i - F_p where A is
A = beta_pos.*(1+ alpha_neg./alpha_pos) ./ (beta_neg - beta_pos);
K_i22 = A.*K_i - alpha_pos - alpha_neg;

%And further K_i2 = F_p*(B - 1) where B is
B = (alpha_pos.*beta_pos + alpha_neg.*beta_pos) ./ ...
    (alpha_neg.*beta_neg + alpha_pos.*beta_pos);
K_i222 =  F_p.*(B - 1);

figure;
subplot(1,2,1);
plot(K_i2, K_i22, 'r.'); axis equal;
xlabel('K_i^*');
ylabel('K_i^{**}');
title('Derived K_i when we set \alpha,\beta_{+/-} = \alpha,\beta_{-/+}')
axis([0 100 0 100]);

subplot(1,2,2);
plot(K_i2, K_i222, 'r.'); axis equal;
xlabel('K_i^*');
ylabel('K_i^{***}');
title('Derived K_i when we set \alpha,\beta_{+/-} = \alpha,\beta_{-/+}')
axis([0 100 0 100]);
%%
%From this we can see K_i2 < 0 iff B < 1, and thus K_i2 < 0 iff
% a_+.b_+ + a_-.b_+ < a_+.b_+ + a_-.b_-
%
% Cancelling out a_+.b_+ from either side of the inequality we have
% a_-.b_+ < + a_-.b_-
%
% Which, given a_- > 0, is true iff b_+ < b_-
%%
%Go through our existing data and check how often Kpos > Kneg in actual
%fitted data. Thankfully very rarely!
for i_sub = 1:10
    study_dir = 'C:\isbe\qbi\data\milano_primovist\';
    visit_dir = [study_dir 'PRIMDCE_' num2str(i_sub) ...
        '\visit1\mdm_analysis24i_reg_DIIRF\'];
    Fpos = load_img_volume([visit_dir 'Fpos.hdr']);
    Fneg = load_img_volume([visit_dir 'Fneg.hdr']);
    Kpos = load_img_volume([visit_dir 'Kpos.hdr']);
    Kneg = load_img_volume([visit_dir 'Kneg.hdr']);
    %
    active_params_model_to_phys(Fpos, Fneg, Kpos, Kneg);
    %
    fitted_voxels = Fneg > 0;
    bad_voxels = Kpos(fitted_voxels) > Kneg(fitted_voxels); 
    num_voxels = sum(fitted_voxels(:));
    num_bad = sum(bad_voxels);
    fprintf('Pt%d: %d bad voxels out of %d\n', i_sub, num_bad, num_voxels);
end


    