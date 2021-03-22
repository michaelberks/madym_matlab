function [M0, T1, Rsq, M0_seed, T1_seed] = fit_spgr_eqn(flipangles, signals, ...
    TRms, bounds, plotoption, fittingmethod)
%FIT_SPGR_EQN returns T1 fit to a spoiled gradient echo sequence
%  [M0, T1] = fit_spgr_eqn(flipangles, signals, TR) returns the T1 (s) and M0
%  value which provides the best fit to the SPGR equation given the
%  variable flip angles (flipangles, degrees) and the associated signals
%  (signals, a.u.) for a repetition time of TR (ms).
%  [M0, T1] = fit_spgr_eqn(flipangles, signals, TR, bounds) limits the T1 fit
%  to bounds(1)<T1<bounds(2).
%  [M0, T1] = fit_spgr_eqn(flipangles, signals, TR, bounds, plotoption) - passing
%  'plot' gives a plot of the fit
%  [M0, T1] = fit_spgr_eqn(flipangles, signals, ... fitmethod) - passing
%  'simplex' or 'lsqcurvefit' allows choice of fitting method
%  A linear approximation is used to seed the fit.  This is followed by a
%  nonlinear estimation to return the final result.
%
%  RAL 5/11/18
%  DJM 19/11/18 edits - return Rsq of fit
%  DJM 30/01/19 edits - add option to plot fit and use simplex method

% Convert TR to seconds
TR = TRms/1000;

% Fitting method -  default to lsqcurevfit if not specified
if exist('fittingmethod','var')
    fitmethod=fittingmethod;
else
    fitmethod='lsqcurvefit';
end

% Fit vfa using a straight line estimate to find an overall T1 and M0
y = signals./sind(flipangles);
x = signals./tand(flipangles);
% fit to the straight line to find p = [m c];
p = polyfit(x,y,1); m = p(1); c = p(2);
M0_seed = abs(c./(1-m));
T1_seed = abs(-TR./log(m));

% Organise bounds
if exist('bounds', 'var')
    lobounds = [0 bounds(1)];
    hibounds = [Inf bounds(2)];
    if T1_seed<bounds(1), T1_seed = bounds(1); end
    if T1_seed>bounds(2), T1_seed = bounds(2); end
else
    lobounds = [0 0];
    hibounds = [Inf Inf];
end

% SPGR equations
spgrEq = @(M0,T1,VFA) M0.*sind(VFA).*(1-exp(-TR./T1))./(1-exp(-TR./T1).*cosd(VFA));
spgr_for_fitting = @(p_nonlin, VFA) spgrEq(p_nonlin(1), p_nonlin(2), VFA);

switch fitmethod
    case 'lsqcurvefit'
        p_nonlin = lsqcurvefit(spgr_for_fitting, [M0_seed, T1_seed], ...
            flipangles, signals, lobounds, hibounds, ...
            optimoptions('lsqcurvefit','Display', 'off'));
        M0 = p_nonlin(1);
        T1 = p_nonlin(2); % we have calculated using TR (s) so T1 is in seconds
    case 'simplex'
        options=optimset('Display', 'off');
        p_nonlin = fminsearchbnd(@objfn,[M0_seed, T1_seed],...
            lobounds, hibounds,options);
        M0 = p_nonlin(1);
        T1 = p_nonlin(2); % we have calculated using TR (s) so T1 is in seconds
end

% Output R-squared of fit
fittedSignals = spgrEq(M0,T1,flipangles);
signal_sse = sum((signals(:)-fittedSignals(:)).^2);
ss_diff_from_mean = sum((signals(:) - mean(signals(:))).^2);
Rsq = 1 - (signal_sse./ss_diff_from_mean);

% Plot if needed
if exist('plotoption','var')
    switch plotoption
        case 'plot'
            flipanglesList=linspace(min(flipangles),max(flipangles),1000)';
            fittedSignalsList=spgrEq(M0,T1,flipanglesList);
            figure
            hold on
            plot(flipangles, signals,'kx','markersize',15,'linewidth',2);
            plot(flipanglesList, fittedSignalsList,'k-',...
                'markersize',15,'linewidth',2);
            title(strcat('T1=',num2str(T1.*1e3),'ms, Rsq=',num2str(Rsq)))
            %xlabel(['FA (' char(176),')'])
            xlabel(['FA (degrees)'])
            ylabel('Signal (a.u.)')
            set(gca,'fontsize',18)
        otherwise
    end
end

% Objective function for simplex fitting
    function objfn = objfn(modelparams)
        % Least squares:
        modelSignals=spgrEq(modelparams(1),modelparams(2),flipangles);
        objfn=sum((signals-modelSignals).^2);
    end

end % end fit_spgr_eqn function