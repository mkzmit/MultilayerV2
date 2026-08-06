function figureHandle = plotTwoLayerTGSResults(results,S)
% Each spot figure contains the fit and local sensitivity

    traceCount = numel(results.traces);
    
    % background set to white for visibility
        figureHandle = figure("Visible",S.figureVisible, ...
            "Name",sprintf("TGS spot %d",results.spot));

    layout = tiledlayout(2,traceCount,"TileSpacing","compact","Padding","compact");
    title(layout,sprintf("TGS fit and sensitivity for spot %d",results.spot))

    for index = 1:traceCount
        trace = results.traces(index);

        nexttile(index)
    
        plot(trace.t*1e6,trace.y,"-","LineWidth",1, ...
            "DisplayName","Measured"); hold on
        plot(trace.t*1e6,trace.yfit,"-","LineWidth",1.6, ...
            "DisplayName","Full fit")
        plot(trace.t*1e6,trace.thermalFit,"--","LineWidth",1.2, ...
            "DisplayName","Thermal fit")
        ylabel("Detector signal")
        title(sprintf("%.4g \\mum, %.1f MHz, %d runs", ...
            trace.Lambda*1e6,trace.acousticFrequency/1e6,numel(trace.run)))
        legend("Location","best")
        grid on
        xlabel("Pump-relative time (\mus)")
    end

    nexttile(traceCount+1,[1,traceCount])
    labels = ["alpha_f","alpha_s","R"];
    fitted = logical(results.parameterFitted(:));
    sensitivity = max(results.sensitivity(fitted),realmin);
    bar(categorical(labels(fitted)),sensitivity)
    set(gca,"YScale","log")
    ylabel("Local sensitivity ||dr/dlog_{10}(p)||_2")
    grid on

    set(findall(figureHandle,"-property","FontName"),"FontName","Cambria")
    set(findall(figureHandle,"-property","FontSize"),"FontSize",10)
end
