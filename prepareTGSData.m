function traces = prepareTGSData(S)
% Load and preprocess measurements for model fitting
% Fuction:
%   Reads the grating-spacing calibration file
%   Extracts identifying information from each measurement filename
%   Matches each measurement to its corresponding baseline trace
%   Interpolates a baseline when an exact nominal-spacing baseline is absent
%   Subtracts the POS and NEG baselines independently
%   Forms the differential signal
%   Detects the pump-arrival time and shifts it to t = 0
%   Retains only the fitting interval
%   Averages repeated runs at each nominal spacing and spot
%   Calculates the calibrated grating wavevector
% Inputs:
%   S ~ Settings structure
% Outputs:
%   traces ~ Structure array containing one averaged trace per nominal
%            spacing and spot combination

    calListing = dir(S.calFile);

    if isempty(calListing)
        error("TGS:Calibration", "No post-processing calibration files were found.");
    end

    runName = strings(0,1);
    LambdaUm = zeros(0,1);

    for i = 1:numel(calListing)
        T = readtable(fullfile(calListing(i).folder, calListing(i).name), ...
            "FileType","text", ...
            "Delimiter"," ", ...
            "MultipleDelimsAsOne",true, ...
            "ReadVariableNames",true, ...
            "VariableNamingRule","preserve");

        runName = [runName; string(T.run_name)];
        LambdaUm = [LambdaUm; double(T.grating_spacing_um)];
    end
    
    % Extract the nominal grating spacing from each calibration run name
        nominalUm = nan(numel(runName),1);
    
        for i = 1:numel(runName)
            beforeSpot = extractBefore(runName(i), "_spot");
            pieces = split(beforeSpot, "_");
            nominalUm(i) = str2double(pieces(end));
        end
    
    % Store calibration 
        cal = containers.Map("KeyType","char", "ValueType","double");
    
        nominalValues = unique(nominalUm(isfinite(nominalUm)));

        for i = 1:numel(nominalValues)
            nominal = nominalValues(i);
            candidates = LambdaUm(nominalUm == nominal & ...
                isfinite(LambdaUm) & LambdaUm > 0);

            if isempty(candidates)
                continue
            end

            [~, closest] = min(abs(candidates - nominal));
            cal(key(nominal)) = candidates(closest)*1e-6;
        end
    
    % Display the final calibration map in ascending nominal order
        calKeys = keys(cal);
        nominalList = str2double(string(calKeys));
        LambdaList = cellfun(@(k) cal(k), calKeys);
        
        [nominalList, order] = sort(nominalList);
        calibratedList = 1e6*LambdaList(order);
        
        disp(table(nominalList(:), calibratedList(:), ...
            'VariableNames', {'Nominal_um','Calibrated_um'}))
    
    
    % Find all files matching the configured data-file pattern
        listing = dir(fullfile(S.dataDir, S.filePattern));
    
    % Preallocate empty structure array containing needed identifiers 
        files = struct( "nominal",{}, "spot",{}, "polarity",{}, "run",{}, ...
            "baseline",{}, "path",{});
    
    for i = 1:numel(listing)
        % Parse the filename using the regular expression supplied in S
            token = regexp(listing(i).name, S.fileRegex, "names", "once");
        % Ignore files that do not follow the expected naming convention
            if isempty(token)
                continue
            end
        % Convert the extracted filename identifiers to their working types
            f.nominal = str2double(token.nominal);
            f.polarity = string(token.polarity);
            f.run = str2double(token.run);
        % A file is a baseline when its parsed location contains "baseline"
            f.baseline = contains(token.location, "baseline");
        
        if f.baseline
            f.spot = NaN;
        else
            f.spot = str2double(extractAfter(string(token.location), "spot"));
        end
    
        f.path = string(fullfile(listing(i).folder, listing(i).name));
    
        % Add the parsed file record to the complete file list.
            files(end+1) = f;
    end
    
    
    
    % Separate ordinary measurement files from baseline
        measurement = find(~[files.baseline]);
        baseline = find([files.baseline]);
    
    pairKey = strings(numel(measurement),1);
    
    for i = 1:numel(measurement)
        f = files(measurement(i));
    
        pairKey(i) = sprintf("%.12g|%d|%d", ...
            f.nominal, f.spot, f.run);
    end
    
    uniquePair = unique(pairKey, "stable");
    
    % Initialize output structure
        traces = struct("Lambda",{}, "q",{}, "t",{}, "y",{}, "tFull",{}, ...
                        "yFull",{}, "t0",{}, "nominal",{}, "spot",{}, ...
                        "run",{}, "files",{});  
    
    for i = 1:numel(uniquePair)
        member = measurement(pairKey == uniquePair(i));
    
        ip = member([files(member).polarity] == "POS");
        in = member([files(member).polarity] == "NEG");

        if numel(ip) ~= 1 || numel(in) ~= 1
            continue
        end
    
        % Use the POS record to obtain the nominal spacing and run identifiers
            f = files(ip);
            k = key(f.nominal);
    
        % Skip the pair if no calibrated spacing is available 
            if ~isKey(cal, k)
                continue
            end
    
        % Find or interpolate the baseline traces for this nominal spacing
            [tbp, ybp, bpPath] = baselineTrace(f.nominal, "POS", ...
                files, baseline, S);
            [tbn, ybn, bnPath] = baselineTrace(f.nominal, "NEG", ...
                files, baseline, S);

        [tp, yp] = readTrace(files(ip).path, S);
        [tn, yn] = readTrace(files(in).path, S);

        % Place each measurement and baseline on a common time grid
            [tp, yp, ybp] = commonGrid(tp, yp, tbp, ybp);
            [tn, yn, ybn] = commonGrid(tn, yn, tbn, ybn);

        % Subtract baseline from each polarity
            yp = yp - ybp;
            yn = yn - ybn;
    
        % Place POS and NEG signals on a common time grid
            [t, yp, yn] = commonGrid(tp, yp, tn, yn);
    
        % Construct differential heterodyne signal
            y = 0.5*(yp - yn);

        t0 = arrivalTime(t, y, S);
        t = t - t0;

        % Retain only fit time interval 
            use = t >= S.tFit(1) & t <= S.tFit(2);
    
        g.Lambda = cal(k);
        g.q = 2*pi/g.Lambda;
        g.t = t(use);
        g.y = y(use);
        g.tFull = [];
        g.yFull = [];
        g.t0 = t0;
        g.nominal = f.nominal;
        g.spot = f.spot;
        g.run = f.run;

        g.files = struct( "POS",files(ip).path, "NEG",files(in).path, ...
            "baselinePOS",bpPath, "baselineNEG",bnPath);
    
        traces(end+1) = g;
    end

    traces = averageRuns(traces,S);
end

function averaged = averageRuns(traces,S)
% Average repeated runs on a shared pump-relative time grid
    groupKey = strings(numel(traces),1);

    for i = 1:numel(traces)
        groupKey(i) = sprintf("%.12g|%d",traces(i).nominal,traces(i).spot);
    end

    uniqueGroup = unique(groupKey,"stable");
    averaged = traces([]);

    for i = 1:numel(uniqueGroup)
        member = traces(groupKey == uniqueGroup(i));
        firstTime = max(arrayfun(@(g) min(g.t),member));
        lastTime = min(arrayfun(@(g) max(g.t),member));
        t = member(1).t;
        t = t(t >= firstTime & t <= lastTime);
        signal = zeros(numel(t),numel(member));

        for j = 1:numel(member)
            signal(:,j) = interp1(member(j).t,member(j).y,t,"linear");
        end

        y = mean(signal,2);
        tFull = t;
        yFull = y;

        % Retain logarithmically distributed points for nonlinear fitting
            if isfield(S,"nFitPoints") && isfinite(S.nFitPoints) && ...
                    S.nFitPoints > 1 && numel(t) > S.nFitPoints
                use = unique(round(logspace(0,log10(numel(t)),S.nFitPoints)));
                t = t(use);
                y = y(use);
            end

        g = member(1);
        g.t = t;
        g.y = y;
        g.tFull = tFull;
        g.yFull = yFull;
        g.t0 = [member.t0];
        g.run = [member.run];
        g.files = [member.files];
        averaged(end+1) = g;
    end
end

function [t,y,path] = baselineTrace(nominal,polarity,files,baseline,S)
% Find an exact baseline or interpolate between bracketing baselines
    records = baseline([files(baseline).polarity] == polarity);
    nominalList = unique([files(records).nominal]);
    exact = records([files(records).nominal] == nominal);

    if isscalar(exact)
        [t,y] = readTrace(files(exact).path,S);
        path = files(exact).path;
        return
    end

    lowerNominal = max(nominalList(nominalList < nominal));
    upperNominal = min(nominalList(nominalList > nominal));

    if isempty(lowerNominal) || isempty(upperNominal)
        error("TGS:Baseline", ...
            "Nominal spacing %.12g is not bracketed by baseline traces.", nominal);
    end

    lower = records([files(records).nominal] == lowerNominal);
    upper = records([files(records).nominal] == upperNominal);

    if numel(lower) ~= 1 || numel(upper) ~= 1
        error("TGS:Baseline", ...
            "Each baseline spacing must contain one %s trace.", polarity);
    end

    [tl,yl] = readTrace(files(lower).path,S);
    [tu,yu] = readTrace(files(upper).path,S);
    [t,yl,yu] = commonGrid(tl,yl,tu,yu);

    weight = (nominal - lowerNominal)/(upperNominal - lowerNominal);
    y = (1 - weight)*yl + weight*yu;
    path = [files(lower).path; files(upper).path];
end

function [t,y] = readTrace(path,S)
% Read and clean the time and signal columns from one data file
    A = readmatrix(path, "NumHeaderLines",S.nHeader);    
    t = A(:,S.tCol)*S.tScale;
    y = A(:,S.yCol);    
    % Retain only rows where both the time and signal are finite
        use = isfinite(t) & isfinite(y);
        t = t(use);
        y = y(use);
end

function [t,yp,yn] = commonGrid(tp,yp,tn,yn)
% Place POS and NEG signals on a shared time grid
    
    if isequal(tp, tn)
        t = tp;
        return
    end    
   
    % Determine the interval over which both polarities contain data
        t1 = max(min(tp), min(tn));
        t2 = min(max(tp), max(tn));        
    
    % Use the time vector with the larger median sample spacing
        if median(diff(tp)) >= median(diff(tn))
            t = tp(tp >= t1 & tp <= t2);
        else
            t = tn(tn >= t1 & tn <= t2);
        end
    
    % Interpolate both polarity signals onto the selected common grid
        yp = interp1(tp, yp, t, "linear");
        yn = interp1(tn, yn, t, "linear");
end

function t0 = arrivalTime(t,y,S)
% Estimate the arrival time of the pump pulse.
    index = find(t >= S.tSearch(1) & t <= S.tSearch(2));
    ys = y(index);
    [~, imax] = max(ys);
    d2y = gradient(gradient(ys));
    [~, id2] = max(d2y(1:imax));
    i0 = max(index(1), index(id2) - S.nShift);
    t0 = t(i0);
end

function k = key(value)
    k = sprintf("%.12g", value);
end
