function fit = fitTwoLayerTGS(traces,S)
% Fit one bounded thermal and acoustic response for all traces in one spot
% The initial thermal fit provides the background used for SAW peak finding
% The final fit combines TwoLayerModel with one damped SAW mode per trace

    fitMask = logical(S.fitParameters);
    xFixed = log10(S.p0);
    x0 = xFixed(fitMask);
    lb = log10(S.pLower(fitMask));
    ub = log10(S.pUpper(fitMask));
    thermalTraces = smoothThermalTraces(traces,S);
    thermalWeights = traceWeights(thermalTraces);
    thermalData = stackedData(thermalTraces,thermalWeights);
    xdata = zeros(size(thermalData));
    options = fittingOptions(S);

    thermalModel = @(x,xdata) stackedThermalModel( ...
        x,thermalTraces,S,thermalWeights);

    xThermal = lsqcurvefit(thermalModel,x0,xdata,thermalData,lb,ub,options);
    [~,thermalState] = stackedThermalModel( ...
        xThermal,thermalTraces,S,ones(size(thermalWeights)));
    [frequency0,tau0] = initialAcousticParameters( ...
        traces,xThermal,thermalState,S);
    fitTraces = finalFitTraces(traces,frequency0,S);

    parameterCount = sum(fitMask);
    traceCount = numel(fitTraces);
    finalX0 = [xThermal,log10(tau0)];
    finalLB = [lb,repmat(log10(S.tauBounds(1)),1,traceCount)];
    finalUB = [ub,repmat(log10(S.tauBounds(2)),1,traceCount)];
    weights = traceWeights(fitTraces);
    ydata = stackedData(fitTraces,weights);
    xdata = zeros(size(ydata));
    model = @(x,xdata) stackedCompositeModel( ...
        x,fitTraces,S,weights,frequency0);

    [x,resnorm,objectiveResidual,exitflag,output,~,J] = ...
        lsqcurvefit(model,finalX0,xdata,ydata,finalLB,finalUB,options);

    [yfit,rebuilt] = stackedCompositeModel( ...
        x,fitTraces,S,ones(size(weights)),frequency0);
    measured = vertcat(fitTraces.y);

    fittedX = xFixed;
    fittedX(fitMask) = x(1:parameterCount);
    fit.x = fittedX;
    fit.fullX = x(:).';
    fit.resnorm = resnorm;
    fit.r = yfit-measured;
    fit.objectiveResidual = objectiveResidual;
    fit.yfit = yfit;
    fit.traces = rebuilt;
    fit.exitflag = exitflag;
    fit.output = output;
    fit.J = J;
    fit.acousticFrequency = frequency0;
    fit.acousticDamping = 10.^x(parameterCount+(1:traceCount));
    fit.parameterFitted = fitMask;
    fit.atLowerBound = false(1,3);
    fit.atUpperBound = false(1,3);
    fit.atLowerBound(fitMask) = ...
        abs(fit.x(fitMask)-lb) <= S.dx;
    fit.atUpperBound(fitMask) = ...
        abs(fit.x(fitMask)-ub) <= S.dx;

    degreesOfFreedom = max(numel(objectiveResidual)-numel(x),1);
    covariance = (resnorm/degreesOfFreedom)*pinv(full(J.'*J));
    xError = sqrt(max(diag(covariance),0));
    p = 10.^fit.x;
    fit.parameterError = nan(1,3);
    fit.parameterError(fitMask) = log(10)*p(fitMask).* ...
        xError(1:parameterCount).';

    % J(:,j)=dr/dx_j
        fit.sensitivity = nan(1,3);
        fit.sensitivity(fitMask) = vecnorm(J(:,1:parameterCount),2,1);
end

function fitTraces = finalFitTraces(traces,frequency,S)
% Begin the final fit after the selected acoustic null-point interval
    fitTraces = traces;

    for g = 1:numel(traces)
        [~,maximum] = max(traces(g).y);
        startTime = traces(g).t(maximum) + ...
            S.fitStartHalfPeriods/(2*frequency(g));
        use = traces(g).t >= startTime;

        if sum(use) < 10
            error("TGS:FitWindow","The final fit window is too short.");
        end

        fitTraces(g).t = traces(g).t(use);
        fitTraces(g).y = traces(g).y(use);
    end
end

function options = fittingOptions(S)
% Create the bounded nonlinear least-squares settings used by both stages
    options = optimoptions("lsqcurvefit", ...
        "Algorithm","trust-region-reflective", ...
        "FiniteDifferenceType","forward", ...
        "FiniteDifferenceStepSize",S.dx, ...
        "MaxIterations",S.maxIter, ...
        "MaxFunctionEvaluations",S.maxEval, ...
        "FunctionTolerance",S.ftol, ...
        "StepTolerance",S.xtol, ...
        "OptimalityTolerance",S.gtol, ...
        "Display",S.display);
end

function thermalTraces = smoothThermalTraces(traces,S)
% Build the smooth background traces used by the initial thermal fit
    thermalTraces = traces;

    for g = 1:numel(traces)
        window = max(3,round(S.smoothTime/median(diff(traces(g).t))));
        thermalTraces(g).y = movmean(traces(g).y,window);
    end
end

function [frequency0,tau0] = initialAcousticParameters( ...
        traces,xThermal,thermalState,S)
% Estimate each SAW frequency from the thermal-fit residual spectrum
    frequency0 = zeros(1,numel(traces));
    tau0 = repmat(S.tau0,1,numel(traces));
    p = physicalParameters(xThermal,S);

    for g = 1:numel(traces)
        t = traces(g).tFull;
        y = traces(g).yFull;
        [Theta,Uz] = TwoLayerModel(t,traces(g).Lambda,p,S);
        D = [Theta,Uz,ones(size(t))];
        background = D*thermalState(g).c;
        residual = y-background;
        dt = median(diff(t));
        derivative = diff(residual)/dt;
        derivative = derivative-mean(derivative);
        nfft = 2^nextpow2(max(2^18,numel(derivative)));
        spectrum = abs(fft(derivative,nfft)).^2;
        frequency = (0:nfft-1).'/(nfft*dt);
        use = frequency >= S.sawFrequencyBounds(1) & ...
            frequency <= S.sawFrequencyBounds(2);
        candidate = frequency(use);
        candidateSpectrum = spectrum(use);
        [~,peak] = max(candidateSpectrum);
        frequency0(g) = candidate(peak);
    end
end

function [yfit,rebuilt] = stackedThermalModel(x,traces,S,weights)
% Stack the fitted thermal model response for all traces in one spot
    p = physicalParameters(x,S);
    n = sum(arrayfun(@(g) numel(g.y),traces));
    yfit = zeros(n,1);
    rebuilt = traces;
    row = 0;

    for g = 1:numel(traces)
        t = traces(g).t;
        y = traces(g).y;
        [Theta,Uz] = TwoLayerModel(t,traces(g).Lambda,p,S);
        D = [Theta,Uz,ones(size(t))];
        [yg,c] = linearFit(D,y);
        rows = row+(1:numel(yg));
        yfit(rows) = weights(g)*yg;
        row = rows(end);
        rebuilt(g).Theta = Theta;
        rebuilt(g).Uz = Uz;
        rebuilt(g).c = c;
        rebuilt(g).yfit = yg;
        rebuilt(g).r = yg-y;
    end
end

function [yfit,rebuilt] = stackedCompositeModel( ...
        x,traces,S,weights,frequency)
% Stack the final thermal plus damped-SAW response for one spot
    traceCount = numel(traces);
    parameterCount = sum(S.fitParameters);
    p = physicalParameters(x(1:parameterCount),S);
    tau = 10.^x(parameterCount+(1:traceCount));
    n = sum(arrayfun(@(g) numel(g.y),traces));
    yfit = zeros(n,1);
    rebuilt = traces;
    row = 0;

    for g = 1:traceCount
        t = traces(g).t;
        y = traces(g).y;
        [Theta,Uz] = TwoLayerModel(t,traces(g).Lambda,p,S);
        decay = exp(-t/tau(g));
        sine = decay.*sin(2*pi*frequency(g)*t);
        cosine = decay.*cos(2*pi*frequency(g)*t);
        D = [Theta,Uz,sine,cosine,ones(size(t))];
        [yg,c] = linearFit(D,y);
        rows = row+(1:numel(yg));
        yfit(rows) = weights(g)*yg;
        row = rows(end);
        rebuilt(g).Theta = Theta;
        rebuilt(g).Uz = Uz;
        rebuilt(g).c = c;
        rebuilt(g).thermalFit = D(:,[1,2,5])*c([1,2,5]);
        rebuilt(g).acousticFit = D(:,3:4)*c(3:4);
        rebuilt(g).yfit = yg;
        rebuilt(g).r = yg-y;
        rebuilt(g).acousticFrequency = frequency(g);
        rebuilt(g).acousticDamping = tau(g);
    end
end

function p = physicalParameters(x,S)
% Combine fitted parameters with values held fixed at their supplied inputs
    fullX = log10(S.p0);
    fullX(logical(S.fitParameters)) = x;
    p = 10.^fullX;
end

function [yfit,c] = linearFit(D,y)
% Solve the linear amplitudes after scaling each model component
    scale = vecnorm(D,2,1);

    if any(~isfinite(D),"all") || any(~isfinite(scale)) || any(scale == 0)
        error("TGS:Model","The bounded model response is not finite.");
    end

    Ds = D./scale;
    cs = Ds\y;
    c = cs./scale.';
    yfit = D*c;
end

function weights = traceWeights(traces)
% Give every grating spacing equal relative weight in the fit
    weights = zeros(numel(traces),1);

    for g = 1:numel(traces)
        scale = norm(traces(g).y-mean(traces(g).y));

        if ~isfinite(scale) || scale == 0
            error("TGS:DataScale","Trace %d has no finite signal variation.",g);
        end

        weights(g) = 1/scale;
    end
end

function ydata = stackedData(traces,weights)
% Stack the weighted measured data for MATLAB nonlinear least squares
    n = sum(arrayfun(@(g) numel(g.y),traces));
    ydata = zeros(n,1);
    row = 0;

    for g = 1:numel(traces)
        rows = row+(1:numel(traces(g).y));
        ydata(rows) = weights(g)*traces(g).y;
        row = rows(end);
    end
end
