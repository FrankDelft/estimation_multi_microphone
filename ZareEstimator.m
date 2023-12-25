    % estimatedSignal: The estimated signal using the weighted average

    [numSamples, numMics] = size(Data);
    weights = ones(numMics, 1) / numMics; % Equal weights
    estimatedSignal = Data * weights; % Weighted average across microphones
    %plot(Data)
    plot(estimatedSignal)

%% Variance calculations

    % Parameters for framing
    frameLength = 0.02; % 20 ms frames
    frameShift = 0.01; % 50% overlap (10 ms)
    frameSize = floor(frameLength * fs);
    shiftSize = floor(frameShift * fs);

    % Number of frames
    numFrames = floor((length(Clean) - frameSize) / shiftSize) + 1;

    % Initialize variables
    K = frameSize; % Frequency bins
    L = numFrames; % Number of time frames
    varemp = zeros(nrmics, 1);

    % Loop over number of microphones
    for m = 1:nrmics
        varianceSum = 0;

        % Process each frame
        for l = 1:L
            frameStart = (l-1) * shiftSize + 1;
            frameEnd = frameStart + frameSize - 1;

            % Frame from clean signal
            cleanFrame = Clean(frameStart:frameEnd);

            % Frame from each microphone and averaging
            micFrame = mean(Data(frameStart:frameEnd, 1:m), 2);

            % FFT of frames
            S = fft(cleanFrame, K);
            EstimatedS = fft(micFrame, K);

            % Variance for this frame
            varianceSum = varianceSum + sum(abs(EstimatedS - S).^2);
        end

        % Average variance for this number of microphones
        varemp(m) = varianceSum / (K * L);
    end

    
