function signal = filter_fir(signal,fs,passType,cutoff,transBand)

nFrames = 1000;

% filter order
m  = firwsord('hamming', fs, transBand);  

% filter coefficients
b  = firws(m, cutoff / (fs / 2), passType, windows('hamming', m + 1)); 

groupDelay = (length(b) - 1) / 2;
dcArray = [1 length(signal) + 1];

for iDc = 1:length(dcArray)-1

    % Pad beginning of data with DC constant and get initial conditions
    ziDataDur = min(groupDelay, dcArray(iDc + 1) - dcArray(iDc));
    [~, zi] = filter(b, 1, double([signal(ones(1, groupDelay) * dcArray(iDc)) ...
        signal(dcArray(iDc):(dcArray(iDc) + ziDataDur - 1))]), [], 2);

    blockArray = [(dcArray(iDc) + groupDelay):nFrames:(dcArray(iDc + 1) - 1) dcArray(iDc + 1)];
    for iBlock = 1:(length(blockArray) - 1)

        % Filter the data
        [signal(1, (blockArray(iBlock) - groupDelay):(blockArray(iBlock + 1) - groupDelay - 1)), zi] = ...
            filter(b, 1, double(signal(1, blockArray(iBlock):(blockArray(iBlock + 1) - 1))), zi, 2);
    end

    % Pad end of data with DC constant
    temp = filter(b, 1, double(signal(1, ones(1, groupDelay) * (dcArray(iDc + 1) - 1))), zi, 2);
    signal(1, (dcArray(iDc + 1) - ziDataDur):(dcArray(iDc + 1) - 1)) = temp(:, (end - ziDataDur + 1):end);

end

