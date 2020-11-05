%So when I refactored the original Madym code, I noticed that compared to
%what seemed like the "correct" implementation of discrete convolution,
%using the factoring trick with exponentials to compute everything in a
%single forward loop, there was a rogue e_delta term used in the average
%between the current and previous points in the AIF. My hunch was this was
%in there as it provides a better "trapezium rule" approximateion to the
%convolution in data with limited temporal resolution. This is confirmed
%with the numerical tests below.
%%
for n_t = [360000 90 360 3600 36000 ]
    
    t = linspace(0,6,n_t);
    if n_t == 360000
        t_full = t;
        Ca_full = population_aif(t_full, 40000);
        Ca_t = Ca_full;
    else
        Ca_t = interp1(t_full, Ca_full, t, 'linear');
    end
    
    k_ep = 100;
    et = exp(-t*k_ep);
    et(isnan(et))=0;
    
    delta_t = t(2)-t(1);
    C_t = zeros(n_t, 4);
    Ct_c = conv(et, Ca_t);
    C_t(:,4) = Ct_c(1:n_t) .* delta_t;
    if n_t == 360000    
        Ct_full = C_t(:,4);
    end
    %
    Ca_t0 = 0;
    e_i = 0;
    for i_t = 2:n_t

        %Get current time, and time change
        t1 = t(i_t);
        delta_t = t1 - t(i_t-1);

        %Compute (offset) combined arterial and vascular input for this time
        Ca_ti = Ca_t(i_t);

        %Update the exponentials for the transfer terms in the two compartments
        e_delta = exp(-delta_t * k_ep);

        %Combine the two compartments with the rate constant to get the final
        %concentration at this time point
        A(1) = delta_t * 0.5 * (Ca_ti + Ca_t0.*e_delta);
        A(2) = delta_t * 0.5 * (Ca_ti + Ca_t0);
        A(3) = delta_t *Ca_ti;
        for i_a = 1:3
            C_t(i_t,i_a) = C_t(i_t-1,i_a).* e_delta + A(i_a);
        end
        Ca_t0 = Ca_ti;
    end
    
    %
    figure; hold all;
    plot(t_full, Ct_full, 'k', 'linewidth', 4.0)
    plot(t, C_t(:,1), 'b');
    plot(t, C_t(:,2), 'r--');
    plot(t, C_t(:,3), 'k-');
    plot(t, C_t(:,4), 'g-.');
    %
    C_ti = interp1(t_full, Ct_full, t, 'linear')';
    display(mean(abs(C_t - C_ti)))
end
    