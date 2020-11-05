%This script checks that our functions for converting between the model
%parameters fitted in the 2CXM and the physiological parameters from which
%they're derived are self-consistent

%% Check phys -> model -> phys
%Generate random values for F_p, PS, v_e and v_p
F_p0 = rand(3,3);
PS0 = rand(3,3);
v_e0 = rand(3,3);
v_p0 = rand(3,3);

%Convert to K_pos, K_neg, F_pos and F_neg
[K_pos, K_neg, F_pos, F_neg] = ...
    two_cx_params_phys_to_model(F_p0, PS0, v_e0, v_p0);

%Convert back to physiological params
[F_p1, PS1, v_e1, v_p1] = ...
    two_cx_params_model_to_phys(K_pos, K_neg, F_pos, F_neg);

%Check differences and confirm all tiny
fprintf('***********************************************************\n');
fprintf('************** Check phys -> model -> phys ****************\n');
fprintf('***********************************************************\n');
fprintf('Change in Fp\n'); 
display(F_p1 - F_p0);

fprintf('Change in PS\n'); 
display(PS1 - PS0);

fprintf('Change in Ve\n'); 
display(v_e1 - v_e0);

fprintf('Change in Vp\n'); 
display(v_p1 - v_p0);

delta_F_p = max(abs(F_p1(:) - F_p0(:)));
delta_PS = max(abs(PS1(:) - PS0(:)));
delta_v_e = max(abs(v_e1(:) - v_e0(:)));
delta_v_p = max(abs(v_p1(:) - v_p0(:)));

fprintf('Max differences, Fp: %4.3f PS: %4.3f Ve: %4.3f Vp: %4.3f\n',...
    delta_F_p, delta_PS, delta_v_e, delta_v_p);

if max([delta_F_p, delta_PS, delta_v_e, delta_v_p] < 1e-6)
    fprintf('Success!!\n');
else
    warning('One or more parameters change more than expected tolerance');
end
fprintf('***********************************************************\n');
%% Check model -> phys -> model
%Generate random values for K_pos, K_neg, F_pos and F_neg
K_pos0 = rand(3,3);
K_neg0 = rand(3,3);
F_pos0 = rand(3,3);
F_neg0 = rand(3,3);

%Convert to F_p, PS, v_e and v_p
[F_p, PS, v_e, v_p] = ...
    two_cx_params_model_to_phys(K_pos0, K_neg0, F_pos0, F_neg0);

%Convert back to model params
[K_pos1, K_neg1, F_pos1, F_neg1] = ...
    two_cx_params_phys_to_model(F_p, PS, v_e, v_p);

%Check differences
fprintf('***********************************************************\n');
fprintf('************** Check model -> phys -> model ***************\n');
fprintf('***********************************************************\n');

fprintf('Change in K+\n');
display(K_pos1 - K_pos0);

fprintf('Change in K_\n'); 
display(K_neg1 - K_neg0);

fprintf('Change in F+\n'); 
display(F_pos1 - F_pos0);

fprintf('Change in F_\n'); 
display(F_neg1 - F_neg0);
fprintf('***********************************************************\n');

%Note that we can get flipping above, with pos and neg swapping (this is
%because phys to model is ambiguous)
swap = abs(F_pos1 - F_pos0) > 1e-6;
K_pos2 = K_pos1;
K_neg2 = K_neg1;
F_pos2 = F_pos1;
F_neg2 = F_neg1;

K_pos2(swap) = K_neg1(swap);
K_neg2(swap) = K_pos1(swap);
F_pos2(swap) = F_neg1(swap);
F_neg2(swap) = F_pos1(swap);

%Check differences
fprintf('Change in K+\n');
display(K_pos2 - K_pos0);

fprintf('Change in K_\n'); 
display(K_neg2 - K_neg0);

fprintf('Change in F+\n'); 
display(F_pos2 - F_pos0);

fprintf('Change in F_\n'); 
display(F_neg2 - F_neg0);

delta_K_pos = max(abs(K_pos2(:) - K_pos0(:)));
delta_K_neg = max(abs(K_neg2(:) - K_neg0(:)));
delta_F_pos = max(abs(F_pos2(:) - F_pos0(:)));
delta_F_neg = max(abs(F_neg2(:) - F_neg0(:)));

fprintf('Max differences, K+: %4.3f K_: %4.3f F+: %4.3f F_: %4.3f\n',...
    delta_K_pos, delta_K_neg, delta_F_pos, delta_F_neg);

if max([delta_K_pos, delta_K_neg, delta_F_pos, delta_F_neg] < 1e-6)
    fprintf('Success!!\n');
else
    warning('One or more parameters change more than expected tolerance');
end
fprintf('***********************************************************\n');
%% Check phys -> model -> phys using F_p and E_neg instead of F_pos, F_neg
%Generate random values for F_p, PS, v_e and v_p
F_p0 = rand(3,3);
PS0 = rand(3,3);
v_e0 = rand(3,3);
v_p0 = rand(3,3);

%Convert to K_pos, K_neg, F_pos and F_neg
[K_pos, K_neg, F_pos, F_neg] = ...
    two_cx_params_phys_to_model(F_p0, PS0, v_e0, v_p0, true);

%Convert back to physiological params
[F_p1, PS1, v_e1, v_p1] = ...
    two_cx_params_model_to_phys(K_pos, K_neg, F_pos, F_neg, true);

%Check differences and confirm all tiny
fprintf('***********************************************************\n');
fprintf('************** Check phys -> model -> phys ****************\n');
fprintf('***********************************************************\n');
fprintf('Change in Fp\n'); 
display(F_p1 - F_p0);

fprintf('Change in PS\n'); 
display(PS1 - PS0);

fprintf('Change in Ve\n'); 
display(v_e1 - v_e0);

fprintf('Change in Vp\n'); 
display(v_p1 - v_p0);

delta_F_p = max(abs(F_p1(:) - F_p0(:)));
delta_PS = max(abs(PS1(:) - PS0(:)));
delta_v_e = max(abs(v_e1(:) - v_e0(:)));
delta_v_p = max(abs(v_p1(:) - v_p0(:)));

fprintf('Max differences, Fp: %4.3f PS: %4.3f Ve: %4.3f Vp: %4.3f\n',...
    delta_F_p, delta_PS, delta_v_e, delta_v_p);

if max([delta_F_p, delta_PS, delta_v_e, delta_v_p] < 1e-6)
    fprintf('Success!!\n');
else
    warning('One or more parameters change more than expected tolerance');
end
fprintf('***********************************************************\n');

%% Check model -> phys -> model using F_p and E_neg instead of F_pos, F_neg
%Generate random values for K_pos, K_neg, F_pos and F_neg
K_pos0 = rand(3,3);
K_neg0 = rand(3,3);
F_p0 = rand(3,3);
E_neg0 = rand(3,3);

%Convert to F_p, PS, v_e and v_p
[F_p, PS, v_e, v_p] = ...
    two_cx_params_model_to_phys(K_pos0, K_neg0, F_p0, E_neg0, true);

%Convert back to model params
[K_pos1, K_neg1, F_p1, E_neg1] = ...
    two_cx_params_phys_to_model(F_p, PS, v_e, v_p, true);

%Check differences
fprintf('***********************************************************\n');
fprintf('************** Check model -> phys -> model ***************\n');
fprintf('***********************************************************\n');

fprintf('Change in K+\n');
display(K_pos1 - K_pos0);

fprintf('Change in K_\n'); 
display(K_neg1 - K_neg0);

fprintf('Change in Fp\n'); 
display(F_p1 - F_p0);

fprintf('Change in E_\n'); 
display(E_neg1 - E_neg0);
fprintf('***********************************************************\n');

%Note that we can get flipping above, with pos and neg swapping (this is
%because phys to model is ambiguous)
swap = abs(K_pos1 - K_pos0) > 1e-6;
K_pos2 = K_pos1;
K_neg2 = K_neg1;
F_p2 = F_p1;
E_neg2 = E_neg1;

K_pos2(swap) = K_neg1(swap);
K_neg2(swap) = K_pos1(swap);
E_neg2(swap) = 1 - E_neg1(swap);

%Check differences
fprintf('Change in K+\n');
display(K_pos2 - K_pos0);

fprintf('Change in K_\n'); 
display(K_neg2 - K_neg0);

fprintf('Change in Fp\n'); 
display(F_p2 - F_p0);

fprintf('Change in E_\n'); 
display(E_neg2 - E_neg0);

delta_K_pos = max(abs(K_pos2(:) - K_pos0(:)));
delta_K_neg = max(abs(K_neg2(:) - K_neg0(:)));
delta_F_p = max(abs(F_p2(:) - F_p0(:)));
delta_E_neg = max(abs(E_neg2(:) - E_neg0(:)));

fprintf('Max differences, K+: %4.3f K_: %4.3f F+: %4.3f F_: %4.3f\n',...
    delta_K_pos, delta_K_neg, delta_F_p, delta_E_neg);

if max([delta_K_pos, delta_K_neg, delta_F_p, delta_E_neg] < 1e-6)
    fprintf('Success!!\n');
else
    warning('One or more parameters change more than expected tolerance');
end
fprintf('***********************************************************\n');

    