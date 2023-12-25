    % estimatedSignal: The estimated signal using the weighted average

    [numSamples, numMics] = size(Data);
    weights = ones(numMics, 1) / numMics; % Equal weights
    estimatedSignal = Data * weights; % Weighted average across microphones
    %plot(Data)
    plot(estimatedSignal)

    
