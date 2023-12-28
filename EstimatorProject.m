%
clear all
close all
clc
% Estimation and Detection Assignment Script
%% Load data
load('Data.mat'); % Assuming Data.mat contains 'Data' and 'Clean'
%%
% Parameters
fs = 16000; % Sampling frequency 
varempEst = zeros(nrmics, 1); % first value of varempEst is very high so the plot will be ugly because of that 
for nrmicsManual = 1:nrmics
%%
frameLength = 0.02; % 20 ms frame
frameShift = 0.01; % 50% overlap (10 ms)

% Convert frame length and overlap to samples
frameSize = floor(frameLength * fs);
shiftSize = floor(frameShift * fs);
stepSize = frameSize - shiftSize;

% Number of frames
numFrames = floor((length(Clean) - frameSize) / shiftSize) + 1;


% Create Hann window
K = frameSize;
L = numFrames; % Number of time frames
window = hann(K);

% Initialize FFT result matrix
fftResult = zeros(frameSize, numFrames,nrmicsManual);
stackedS = zeros(frameSize, numFrames);
for m = 1:nrmicsManual
    % Divide signal into frames and perform FFT
    for i = 1:numFrames
        frameStart = (i-1) * stepSize + 1;
        frameEnd = frameStart + frameSize - 1;
        % Apply Hann window and FFT
        windowedFrame = Data(frameStart:frameEnd,m) .* window;
        fftResult(:, i, m) = fft(windowedFrame);
    end
end
   
%%
%estimate the noise in the first second

%now to make an estimate of the noise based on the 1st second of noise
Noise_duration=1;
L_1s=(Noise_duration-frameLength)/(frameLength*(0.5))+1;
%variance over L frames
estimated_noise_psd=zeros(nrmicsManual,K);

for m = 1:nrmicsManual
    mic_m_sum=zeros(K,1);
    for l = 1:L_1s
        frameStart = (l-1) * shiftSize + 1;
        frameEnd = frameStart + frameSize - 1;
        
        noisyFrame=Data(frameStart:frameEnd,m).*window;
        noise_FFT=fft(noisyFrame,K);
        mic_m_sum=mic_m_sum+abs(noise_FFT).^2;
        
    end
    estimated_noise_psd(m,:)=(mic_m_sum').*1/L_1s;
end
%% average over the frequency bins and 
var_m=mean(estimated_noise_psd,2);


%% 

% Estimate PSDs for signal and noise
Sxx = zeros(K, L);
Nxx = zeros(K, L);
% for m = 1:nrmics
%     for i = 1:numFrames
%         frameStart = (i-1) * stepSize + 1;
%         frameEnd = frameStart + frameSize - 1;
% 
%         % Clean and noisy frames
%         cleanFrame = Clean(frameStart:frameEnd) .* window;
%         noisyFrame = Data(frameStart:frameEnd, m) .* window;
% 
%         % FFT
%         cleanPSD = abs(fft(cleanFrame)).^2;
%         noisyPSD = abs(fft(noisyFrame)).^2;
% 
%         % PSD estimation
%         Sxx(:, i) = Sxx(:, i) + cleanPSD / nrmics;
%         noisePSD = noisyPSD - cleanPSD;
%         Nxx(:, i) = Nxx(:, i) + noisePSD / nrmics;
%     end
% end


for m = 1:nrmicsManual
    for i = 1:numFrames
        frameStart = (i-1) * stepSize + 1;
        frameEnd = frameStart + frameSize - 1;

        % Clean and noisy frames
        noisyFrame = Data(frameStart:frameEnd, m) .* window;

        % FFT
        cleanPSDest = abs(fft(noisyFrame)).^2-estimated_noise_psd(m,:)';
        noisePSD = estimated_noise_psd(m,:)';

        % PSD estimation
        Sxx(:, i) = Sxx(:, i) + cleanPSDest / nrmicsManual;
        Nxx(:, i) = Nxx(:, i) + noisePSD / nrmicsManual;
    end
end

%%
% Initialize W
W = complex(zeros(K, L, nrmicsManual));

% Calculate Wiener filter weights
for k = 1:K
    for l = 1:L
        % Wiener filter
        W(k, l, :) = Sxx(k, l) ./ (Sxx(k, l) + Nxx(k, l));
    end
end

%%


% Initialize W with complex values
%W = complex(rand(frameSize, numFrames, nrmics), rand(frameSize, numFrames, nrmics));

% Perform the multiplication and sum
for k = 1:K
    for l = 1:L
        % Element-wise multiplication and sum across the M dimension
        stackedS(k, l) = sum(conj(W(k, l, :)) .* fftResult(k, l, :))/nrmicsManual;
    end
end
  varianceSum = 0;

        % Process each frame
        for l = 1:L
            frameStart = (l-1) * shiftSize + 1;
            frameEnd = frameStart + frameSize - 1;

            % Frame from clean signal
            cleanFrame = Clean(frameStart:frameEnd);

            % FFT of frames
            S = fft(cleanFrame.*window, K);

            % Variance for this frame
            varianceSum = varianceSum + sum(abs(stackedS(:,l) - S).^2);
        end

        % Average variance for this number of microphones
        varempEst(nrmicsManual) = varianceSum / (K * L);
end
%%
%now lets get the fischer information for different numbers of microphones
cum_sum_inverse_estimated_noise=cumsum(1./var_m);

%now lets calculate the CRLB for different numbers of microphones
crlb_m=1./cum_sum_inverse_estimated_noise;

plot(crlb_m)
hold on;
plot(varempEst)
hold off;
%legend("crlb","actual")
   


%%
reconstructedSignal = zeros((numFrames-1) * stepSize + frameSize, 1);

% Overlap-add method
for i = 1:numFrames
    frameStart = (i-1) * stepSize + 1;
    frameEnd = frameStart + frameSize - 1;

    % Inverse FFT
    ifftFrame = ifft(stackedS(:, i), 'symmetric');

    % Overlap-add
    reconstructedSignal(frameStart:frameEnd) = reconstructedSignal(frameStart:frameEnd) + ifftFrame;
end
audiowrite("reconstructed.wav", reconstructedSignal, 16000);
%%
figure;
plot(Clean, 'b', 'DisplayName', 'Original');
hold on;
plot(reconstructedSignal, 'r', 'DisplayName', 'Modeled');


xlabel('Time');
ylabel('Value');
title('Original vs Modeled Time Series');
legend('show');
grid on;
