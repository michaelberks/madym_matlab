%% AUEM
F_p = 0.5;
v_ecs = 0.15;
k_i = 0.07;
k_ef = 0.05;

[Apos, Aneg, Bpos, Bneg] = active_params_phys_to_model(...
    F_p, v_ecs, k_i, k_ef);
[Apos3, Aneg3, Bpos3, Bneg3] = active_params_phys_to_model(...
    F_p, v_ecs, k_i, 0);
[Apos2, Aneg2, Bpos2, Bneg2] = active_params_phys_to_model(...
    F_p, v_ecs, 0, 0);

display([Apos, Aneg, Bpos, Bneg]);
display([Apos3, Aneg3, Bpos3, Bneg3]);
display([Apos2, Aneg2, Bpos2, Bneg2]);
%%
F_p2 = Apos + Aneg;
display([F_p F_p2]);
%%
Ei = k_i / (k_i + F_p);
Epos = Apos / (Apos + Aneg);
%%
B = 1 - Bpos/Bneg;
Ei2 = Epos*B;
display([Ei Ei2]);
%%
k_i2 = Ei*F_p / (1 - Ei);
display([k_i k_i2]);
%%
k_i3 = Apos*B / (1 - Ei);
display([k_i k_i3]);
%%
Ein = 1 - Ei;
Ein2 = 1 - Epos*B;
Ein3 = (Apos*(1-B) + Aneg) / (Apos + Aneg);
display([Ein Ein2 Ein3]);
%%
k_i4 = (Apos + Aneg)*Apos*B / (Apos*(1-B) + Aneg);
display([k_i k_i4]);
%%
k_i5 = (Apos + Aneg)*Apos*(Bneg-Bpos) / (Apos*Bpos + Aneg*Bneg);
display([k_i k_i5]);
%%
v_ecs2 = (F_p + k_i) / Bneg;
display([v_ecs v_ecs2]);
%%
v_ecs3 = (Apos+Aneg)/Bneg + (Apos + Aneg)*Apos*(Bneg-Bpos) / (Bneg*(Apos*Bpos + Aneg*Bneg));
display([v_ecs v_ecs3]);
%%
v_ecs4 = (Apos+Aneg)/Bneg * (1 + Apos*(Bneg-Bpos)/(Apos*Bpos + Aneg*Bneg));
display([v_ecs v_ecs4]);
%%
v_ecs5 = (Apos + Aneg)^2 / (Apos*Bpos + Aneg*Bneg);
display([v_ecs v_ecs5]);
%%
k_ef2 = (1 - v_ecs)*Bpos;
display([k_ef k_ef2]);
%%
k_ef3 = Bpos - Bpos*(Apos + Aneg)^2 / (Apos*Bpos + Aneg*Bneg);
display([k_ef k_ef3]);
%% ------------------------------------------------------------------------
%% 2CXM
%% ------------------------------------------------------------------------
F_p = 0.6;
PS = 0.1;
v_e = 0.2;
v_p = 0.05;
[Apos, Aneg, Bpos, Bneg] = two_cx_params_phys_to_model(...
    F_p, PS, v_e, v_p);
display([Apos, Aneg, Bpos, Bneg]);
%%
% From previous work, v_p = v_ecs
v_p2 = (Apos + Aneg)^2 / (Apos*Bpos + Aneg*Bneg);
display([v_p v_p2]);
%%
v_e2 = (Apos+Aneg)/Bpos + Aneg*(Bpos-Bneg)/(Bpos.*Bneg)...
    - (Apos+Aneg)^2 / (Bpos*Apos + Bneg*Aneg);
display([v_e v_e2]);
%%
v_e3 = (Aneg*Bpos + Apos*Bneg) / (Bpos*Bneg) - (Apos+Aneg).^2 ./ (Bpos.*Apos + Bneg.*Aneg);
display([v_e v_e3]);

v_e4 = Apos*Aneg*(Bpos-Bneg)^2 / (Bpos*Bneg*(Apos*Bpos + Aneg*Bneg));
display([v_e v_e4]);
%%
Epos = Apos / F_p;
PS2 = F_p/(Epos.*(Bpos-Bneg)+Bneg)  * (Bpos - Bpos*Bneg /(Epos*(Bpos-Bneg)+Bneg) - Epos*(Bpos-Bneg));
display([PS PS2]);
%%
PS3 = ( (Apos + Aneg)^2 / (Apos*Bpos + Aneg*Bneg) ) *...
    ((Apos*Bneg + Aneg*Bpos) / (Apos+Aneg) - Bpos*Bneg*(Apos+Aneg) / (Apos*Bpos + Aneg*Bneg));


display([PS PS3]);
%%
PS4 = Apos*Aneg*(Apos+Aneg)*(Bpos-Bneg)^2 / (Apos*Bpos + Aneg*Bneg)^2;
display([PS PS3]);