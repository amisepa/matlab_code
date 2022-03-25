%% Matlab displayspectrogram function

function displayspectrogram(t,f,Pxx,isFsnormalized,faxisloc,esttype, threshold)
% Cell array of the standard frequency units strings

if strcmpi(esttype,'power')
    plotOpts.cblbl = getString(message('signal:dspdata:dspdata:PowerdB'));
else
    if isFsnormalized
        plotOpts.cblbl = getString(message('signal:dspdata:dspdata:PowerfrequencydBradsample'));
    else
        plotOpts.cblbl = getString(message('signal:dspdata:dspdata:PowerfrequencydBHz'));
    end
end

%Threshold in dB
plotOpts.freqlocation = faxisloc;
plotOpts.threshold = 10*log10(threshold+eps);
plotOpts.isFsnormalized = logical(isFsnormalized);

%Power in dB
signalwavelet.internal.convenienceplot.plotTFR(t,f, 10*log10(abs(Pxx)+eps),plotOpts);    

