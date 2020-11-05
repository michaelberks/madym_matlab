function [] = active_model_diagram(F_p, v_e, k_i, k_ef, fa, ax)
%ACTIVE_MODEL_DIAGRAM Make a scaled diagram of model contrast flow
%   [] = active_model_diagram(F_p, v_e, k_i, k_ef, fa)
%
% Inputs:
%      F_p, v_e, k_i, k_ef, fa - Inputs to active model, all scalars
%
%      ax - axes to plot on, calls gca if not supplied
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 08-Mar-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('ax', 'var') || isempty(ax)
    gca;
else
    axes(ax);
end
font_size = 14;

%Make boxes for Ve and Vi
v_ex = 100*v_e;
ve_xx = [0 v_ex v_ex 0];
vi_xx = [v_ex 100 100 v_ex];
v_yy = [0 0 100 100];

%Make arterial input arrow
flow_arrow_base_width = 100;
Fa_w = F_p*fa*flow_arrow_base_width;
Fa_xc = v_ex / 2;
Fa_xl = -0.5*Fa_w+Fa_xc;
Fa_xr = 0.5*Fa_w+Fa_xc;
Fa_xx = [Fa_xc Fa_xr Fa_xr Fa_xl Fa_xl];
Fa_yy = [95 105 120 120 105];  

%Make venous input arrow
Fv_h = F_p*(1-fa)*flow_arrow_base_width;
Fv_yc = 70;
Fv_yt = 0.5*Fv_h + Fv_yc;
Fv_yb = -0.5*Fv_h + Fv_yc;
Fv_xx = [-20 -5 5 -5 -20];
Fv_yy = [Fv_yt Fv_yt Fv_yc Fv_yb Fv_yb];

%Make flow output arrow
Fo_w = F_p*flow_arrow_base_width;
Fo_xc = v_ex / 2;
Fo_xl = -0.5*Fo_w+Fo_xc;
Fo_xr = 0.5*Fo_w+Fo_xc;
Fo_xx = [Fo_xc Fo_xr Fo_xr Fo_xl Fo_xl];
Fo_yy = [-20 -10 5 5 -10];

%Make k_i arrow
ki_h = 200*k_i;
ki_yc = 50;
ki_yt = 0.5*ki_h + ki_yc;
ki_yb = -0.5*ki_h + ki_yc;
ki_xx = v_ex+[0 15 25 15 0]-5;
ki_yy = [ki_yt ki_yt ki_yc ki_yb ki_yb];

%Make k_ef arrow
ke_h = 200*k_ef;
ke_yc = 50;
ke_yt = 0.5*ke_h + ke_yc;
ke_yb = -0.5*ke_h + ke_yc;
ke_xx = 95+[0 15 25 15 0];
ke_yy = [ke_yt ke_yt ke_yc ke_yb ke_yb];


axis equal; axis off; hold all;
patch(ve_xx, v_yy, [162 194 194]/255);
patch(vi_xx, v_yy, [194 162 194]/255);
patch(Fa_xx, Fa_yy, 'r');
patch(Fv_xx, Fv_yy, 'b');
patch(Fo_xx, Fo_yy, 'm');
patch(ki_xx, ki_yy, 'c');
patch(ke_xx, ke_yy, 'g');

text(v_ex/2,80, 'Ve', 'fontsize', font_size);
text(50,50, 'Vi', 'fontsize', font_size);
text(Fa_xc, sum(Fa_yy(2:3))/2, 'Fa', 'fontsize', font_size);
text(sum(Fv_xx(1:2))/2, Fv_yc, 'Fv', 'fontsize', font_size);
text(Fo_xc, sum(Fo_yy(2:3))/2, 'Fp', 'fontsize', font_size);
text(v_ex,ki_yc, 'Ki', 'fontsize', font_size);
text(100,ke_yc, 'Kef', 'fontsize', font_size);