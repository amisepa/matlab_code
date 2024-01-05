function EEG = injectHeartArtifacts(EEG, heartRate, artifactAmplitude)
    % EEG - EEGLAB EEG structure
    % heartRate - Heart rate in beats per minute (BPM)
    % artifactAmplitude - Amplitude of the heart artifacts

    % Calculate the interval of heart beats in samples
    heartInterval = 60 / heartRate * EEG.srate; % Convert BPM to sample interval

    % Define a basic shape for the heart artifact (QRS complex)
    qrsWidth = round(0.1 * EEG.srate); % 0.1 second width for QRS complex
    qrsShape = artifactAmplitude * [zeros(1, qrsWidth/2), ones(1, qrsWidth), zeros(1, qrsWidth/2)];

    % Inject artifacts into EEG.data
    for i = 1:size(EEG.data, 3) % Loop through epochs
        for j = 1:size(EEG.data, 1) % Loop through channels
            % Determine the start points for each QRS complex
            qrsStarts = 1:heartInterval:size(EEG.data, 2) - length(qrsShape);
            for startIdx = qrsStarts
                % Add the QRS complex to the EEG data
                EEG.data(j, startIdx:startIdx + length(qrsShape) - 1, i) = ...
                    EEG.data(j, startIdx:startIdx + length(qrsShape) - 1, i) + qrsShape;
            end
        end
    end
end