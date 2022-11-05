function [isFS,isTAN,isSPN,isLowFRThin]=classifyStriatumUnits(halfmaxWidth,meanFR)

% based on Berke, Okatan, Skurski, Eichenbaum. 2004. Neuron.

if halfmaxWidth<2*10^-4
    if meanFR>=10^0.1
        % FS
        isFS=true;
        isTAN=false;
        isSPN=false;
        isLowFRThin=false;
    else
        isFS=false;
        isTAN=false;
        isSPN=false;
        isLowFRThin=true;
    end
else
    if meanFR>4
        % TAN
        isFS=false;
        isTAN=true;
        isSPN=false;
        isLowFRThin=false;
    else
        % SPN
        isFS=false;
        isTAN=false;
        isSPN=true;
        isLowFRThin=false;
    end
end

end