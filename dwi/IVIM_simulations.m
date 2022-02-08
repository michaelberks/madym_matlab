%--------------------------------------------------------------------------
% This script provides a quick demonstration of fitting signals generated 
% for diffusion IVIM model using the Matlab wrapper to Madym's DWI tools
%
% The example shows running a fit at several noise levels for a fixed set
% of IVIM ground truth parameters, however this basic simulation may easily
% be extended to larger sets of more variable inputs, or with different 
% B-values etc.
%%
%Set B-values at which to simulate signal inputs
B_vals = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 300.0, 500.0, 800.0];

%Set ground truth IVIM parameters
S0 = 100;
D = 0.8e-3;
f = 0.2;
D_star = 15e-3;
n_samples = 1e4;

for sigma = [0.1, 1, 2, 4, 10]
    IVIM_simulation(n_samples, B_vals, S0, D, f, D_star, sigma);
end

function IVIM_simulation(n_samples, B_vals, S0, D, f, D_star, sigma)

    %Generate simulated IVIM test data with Rician noise added
    S0s = repmat(S0, n_samples,1);
    Ds = repmat(D, n_samples,1);
    fs = repmat(f,n_samples,1);
    D_stars = repmat(D_star,n_samples,1);

    signals = IVIM_model(B_vals, S0s, Ds, fs, D_stars);
    signals_n = add_rician_noise(signals, sigma);

    %Fit IVIM using Madym
    Bvals_thresh = [40.0,60.0,100.0,150.0];
    [model_params] = run_madym_DWI_lite(...
        'IVIM', signals_n, B_vals,...
        'Bvals_thresh', Bvals_thresh);

    %Plot the output
    param_names = {'S_0', 'D', 'f', 'D^*'};
    param_units = {'a.u.', 'mm^2/s', 'no units', 'mm^2/s'};
    gt = [S0 D f D_star];

    n_bins = round(n_samples / 100);
    figure('Name', sprintf('Rician noise sigma = %d', sigma));
    for i_p = 1:4
        [counts, bins] = hist(model_params(:,i_p), n_bins); %#ok
        med_p = median(model_params(:,i_p));
        max_count = max(counts);
        subplot(2,2,i_p); hold all;
        bar(bins, counts);
        plot([gt(i_p), gt(i_p)], [0 max_count], 'g-', 'linewidth', 2);
        plot([med_p, med_p], [0 max_count], 'r--', 'linewidth', 2);

        legend({...
            'Distribution of fitted parameters',...
            'Ground truth',...
            'Median of fitted parameters'});

        xlabel(sprintf('%s (%s)', param_names{i_p},  param_units{i_p}));
        ylabel(sprintf('Histogram counts per %d', n_samples));
        title(sprintf('Fitting %s = %5.4f: median fit = %5.4f',...
            param_names{i_p}, gt(i_p), med_p));
    end
end
    
    
    
    
    
