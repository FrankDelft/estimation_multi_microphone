    % estimatedSignal: The estimated signal using the weighted average
   
    [numSamples, numMics] = size(Data);
    weights = ones(numMics, 1) / numMics; % Equal weights
    estimatedSignal = Data * weights; % Weighted average across microphones
    
    %plot(Data)
    plot(estimatedSignal)
    fs=16000;
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
    window=hann(K);

    
    
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
            S = fft(cleanFrame.*window, K);
            EstimatedS = fft(micFrame.*window, K);

            % Variance for this frame
            varianceSum = varianceSum + sum(abs(EstimatedS - S).^2);
        end

        % Average variance for this number of microphones
        varemp(m) = varianceSum / (K * L);
    end

    
%% 


%now to make an estimate of the noise based on the 1st second of noise
Noise_duration=1;
L_1s=(Noise_duration-frameLength)/(frameLength*(0.5))+1;
 %variance over L frames
estimated_noise=zeros(nrmics,K);

for m = 1:nrmics
    mic_m_sum=zeros(K,1);
    for l = 1:L_1s
        frameStart = (l-1) * shiftSize + 1;
        frameEnd = frameStart + frameSize - 1;
        
        noise_FFT=fft(Data(frameStart:frameEnd,m).*window,K);
        mic_m_sum=mic_m_sum+abs(noise_FFT).^2;
    end
    estimated_noise(m,:)=(mic_m_sum').*(1/L_1s);
    
end
%average over the frequency bins and 
var_m=mean(estimated_noise,2);

%% 

%now lets get the fischer information for different numbers of microphones
cum_sum_inverse_estimated_noise=cumsum(1./var_m);

%now lets calculate the CRLB for different numbers of microphones
crlb_m=1./cum_sum_inverse_estimated_noise;

plot(crlb_m)
hold on;
plot(varemp)
hold off;
%legend("crlb","actual")




    
