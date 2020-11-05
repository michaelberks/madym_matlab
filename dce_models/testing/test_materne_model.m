%% ------------------------------------------------------------------------
% Script testing the generating fitting of the Materne DCE model
%% ------------------------------------------------------------------------

%% Set default parameters and ranges
param_names = {'k1a', 'k1v', 'k2', 'aoffset', 'voffset', 'all'};
param_titles = {'k_{1a}', 'k_{1v}', 'k_{2}', 't_a', 't_v', 'All'};

params.k1a.def = 0.3;
params.k1v.def = 0.3;
params.k2.def = 0.5;
params.aoffset.def = 0.1;
params.voffset.def = 0.1;

params.k1a.range = linspace(0, 1, 11);
params.k1v.range = linspace(0, 1, 11);
params.k2.range = linspace(0, 5, 11);
params.aoffset.range = linspace(0, 0.5, 11);
params.voffset.range = linspace(0, 0.5, 11);

dyn_times = 5*(0:99)/60;
t_injection = 8;
aif = population_aif(dyn_times, t_injection); 
pif = compute_PIF(aif, [], dyn_times);

T1 = ones(1,11)*600;
S0 = 80;
FA = 18;
TR = 2.4;
relax_coeff = 14e-3;
%% Look at the effect of each parameter in turn
dce_ca = cell(6,1);
for i_param = 1:6
    model_params = zeros(5, 11);
    
    for j_param = 1:5
        if i_param == j_param || i_param == 6
            model_params(j_param,:) = ...
                params.(param_names{j_param}).range;
        else
            model_params(j_param,:) = ...
                ones(1,11)*params.(param_names{j_param}).def;
        end
    end

    dce_ca{i_param} = ...
        materne_model(...
        model_params(1,:), model_params(2,:), model_params(3,:),...
        model_params(4,:), aif, pif, dyn_times, 0, model_params(5,:));
end

%%
figure;
for i_param = 1:6
    subplot(2,3,i_param);
    plot(dyn_times, dce_ca{i_param});
    xlabel('Time (mins)');
    ylabel('C_t');
    title(param_titles{i_param});
end
%%
%Now try and fit to each varying parameter
for i_param = 1:6
    dce_si = concentration_to_signal(...
        dce_ca{i_param}, FA, TR, ...
        T1, S0, relax_coeff, t_injection);
    
    noise = randn(size(dce_si))*2;

    figure; 
    a1 = subplot(1,2,1);
    plot(dyn_times, dce_ca{i_param});
    hold all;
    subplot(1,2,2);
    plot(dyn_times, dce_si + noise);
    
    [fit.k1a, fit.k1v, fit.k2, fit.aoffset, fit.voffset] = fit_materne_model(...
        dce_si+noise, dyn_times, aif, pif, T1,...
        'FA', FA, 'TR', TR, 'num_itr', 100, 'num_baseline', t_injection);
    
    dce_ca_fit = ...
        materne_model(...
        fit.k1a, fit.k1v, fit.k2, fit.aoffset, aif, pif,...
        dyn_times, 0, fit.voffset);
    
    plot(a1, dyn_times, dce_ca_fit, 'k--');
    
    fprintf('***********\n');
    fprintf('Fitting %s\n', param_names{i_param});
    
    for j_param = 1:5
        if i_param == j_param || i_param == 6
            p_i = params.(param_names{j_param}).range;
        else
            p_i = ones(1,11)*params.(param_names{j_param}).def;
        end
        
        fprintf('%s = \n', param_names{j_param});
        fprintf('%f %f\n', [p_i(:)'; fit.(param_names{j_param})(:)']);
    end           
end
%%
dyn_series = zeros(6, 11, 1, 100);
for i_param = 1:6
    dce_si = concentration_to_signal(...
        dce_ca{i_param}, FA, TR, ...
        T1, S0, relax_coeff, t_injection);
    
    dyn_series(i_param, :, :, :) = size(permute(dce_si, [4 2 3 1]));
end
    








