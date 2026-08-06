function results = MainTwoLayerScript()
% Purpose: Fit one shared [alpha_f, alpha_s, R] for each spot
% Instructions: Edit the user input block, then run 

%% USER INPUTS 
    % Files
        S.dataDir = fullfile(fileparts(mfilename("fullpath")), ...
            "MKZ_coatings-20260806T133341Z-1-001", "MKZ_coatings", ...
            "Nb_10klbs", "2026-07-21"); % Path to folder containing data files
        S.calFile = fullfile(S.dataDir, "*_postprocessing.txt"); % Path to post-processing file
        S.filePattern = "*.txt"; 
        S.fileRegex = "(?<nominal>\d+\.\d+)_" + "(?<location>baseline|spot\d+(?:-baseline)?)-" + ...
            "(?<polarity>POS|NEG)-(?<run>\d+)\.txt$"; % change based on naming convention
    % File layout
        S.nHeader = 16; % Number of header lines to skip
        S.tCol = 1; % First column of file is time
        S.yCol = 2; % Second column of file is signal
        S.tScale = 1; % Use if time isn't in seconds
    % Pump arrival and common fit window
        S.tSearch = [0,100e-9]; % Search interval used to locate pump arrival [s]
        S.nShift = 0; % Number of samples placed before the detected arrival
        S.tFit = [0,inf]; % Time interval retained for fitting [s]
        S.nFitPoints = 800; % Maximum number of averaged time points per trace
        S.smoothTime = 4e-9; % Moving-average interval used for the initial thermal fit [s]
    % Known film/substrate quantities
        S.L = 1.5e-6; % film thickness [m]
        S.Cf = 2.56e6 ; % film volumetric heat capacity [J m^-3 K^-1]
        S.Cs = 2.27e6; % substrate volumetric heat capacity [J m^-3 K^-1]
        S.Ef = 4.11e11; % film Young's modulus [Pa]
        S.Es = 1.05e11; % substrate Young's modulus [Pa]
        S.nuf = 0.28; % film Poisson ratio
        S.nus = 0.40; % substrate Poisson ratio
        S.alphathf = 4.42e-6; % film thermal expansion coefficient [K^-1]
        S.alphaths = 7.07e-6; % substrate thermal expansion coefficient [K^-1]
    % Initial nonlinear parameters [alpha_f, alpha_s, R]
        S.p0 = [6.7e-5,2.7e-5,1e-8];  % [m^2/s, m^2/s, m^2 K/W]
        S.pLower = [1.5e-5,1.5e-5,1e-14]; % Lower parameter bounds
        S.pUpper = [1e-1,1e-1,1e-2]; % Upper parameter bounds
        S.fitParameters = [true,true,true]; % Fit [alpha_f, alpha_s, R] for testing use
    % Fourier calculation
        S.Nw = 256; % even number of frequency samples
        S.wmin = 2*pi*1e3; % [rad/s]
        S.wmax = 2*pi*1e12; % [rad/s]
        S.Q0 = 1; % unit absorbed surface impulse
    % Bounded nonlinear least squares
        S.maxIter = 10000;
        S.maxEval = 2000;
        S.ftol = 1e-8;
        S.xtol = 1e-8;
        S.gtol = 1e-8;
        S.dx = 2e-3; % finite-difference step in log10 coordinates
        S.display = "iter";
    % Acoustic fitting
        S.sawFrequencyBounds = [1e8,9e8]; % SAW frequency bounds [Hz]
        S.fitStartHalfPeriods = 2; % Number of acoustic half-periods after the signal maximum
        S.tau0 = 5e-8; % Initial acoustic damping time [s]
        S.tauBounds = [5e-9,7e-7]; % Acoustic damping-time bounds [s]
    % Plotting
        S.makePlots = true;
        S.figureVisible = "on";

%% DATA PREPARATION
    traces = prepareTGSData(S);
    
    if isempty(traces)
        error("TGS:Data","No complete calibrated POS/NEG/baseline sets were found.");
    end

%% FIT EACH SPOT
    spotList = unique([traces.spot]);
    spotResults = cell(numel(spotList),1);

    for i = 1:numel(spotList)
        spot = spotList(i);
        spotTraces = traces([traces.spot] == spot);
        fit = fitTwoLayerTGS(spotTraces,S);

    %% RECONSTRUCTION
        r = fit.r;
        p = 10.^fit.x;
        spotResult.spot = spot;
        spotResult.alpha_f = p(1);
        spotResult.alpha_s = p(2);
        spotResult.R = p(3);
        spotResult.x = fit.x;
        spotResult.sensitivity = fit.sensitivity;
        spotResult.parameterError = fit.parameterError;
        spotResult.traces = fit.traces;
        spotResult.r = r;
        spotResult.resnorm = sum(r.^2);
        spotResult.exitflag = fit.exitflag;
        spotResult.output = fit.output;
        spotResult.atLowerBound = fit.atLowerBound;
        spotResult.atUpperBound = fit.atUpperBound;
        spotResult.parameterFitted = fit.parameterFitted;
        spotResult.acousticFrequency = fit.acousticFrequency;
        spotResult.acousticDamping = fit.acousticDamping;

        spotResult.summary = table( ...
            repmat(spot,3,1),["alpha_f";"alpha_s";"R"],p(:),fit.sensitivity(:), ...
            fit.parameterError(:),fit.parameterFitted(:), ...
            fit.atLowerBound(:),fit.atUpperBound(:), ...
            ["m^2/s";"m^2/s";"m^2 K/W"], ...
            'VariableNames',{'Spot','Symbol','Value','LocalSensitivity', ...
            'StandardError','Fitted','AtLowerBound','AtUpperBound','Unit'});

        fprintf("spot %d\n",spot);
        disp(spotResult.summary)
        fprintf("exitflag = %d\n",fit.exitflag);
        fprintf("iterations = %g\n",fieldValue(fit.output,"iterations"));
        fprintf("function evaluations = %g\n",fieldValue(fit.output,"funcCount"));
        fprintf("||r||_2 = %.6g\n",norm(r));

        if S.makePlots
            spotResult.figure = plotTwoLayerTGSResults(spotResult,S);
        else
            spotResult.figure = gobjects(0);
        end

        spotResults{i} = spotResult;
    end

    results.spots = [spotResults{:}];
    results.summary = vertcat(results.spots.summary);
end


function value = fieldValue(s,name)
if isfield(s,name)
    value = s.(name);
else
    value = NaN;
end
end
