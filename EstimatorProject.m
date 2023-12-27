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
fftResult = zeros(frameSize, numFrames,nrmics);
stackedS = zeros(frameSize, numFrames);
for m = 1:nrmics
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

% Estimate PSDs for signal and noise
Sxx = zeros(K, L);
Nxx = zeros(K, L);
for m = 1:nrmics
    for i = 1:numFrames
        frameStart = (i-1) * stepSize + 1;
        frameEnd = frameStart + frameSize - 1;

        % Clean and noisy frames
        cleanFrame = Clean(frameStart:frameEnd) .* window;
        noisyFrame = Data(frameStart:frameEnd, m) .* window;

        % FFT
        cleanFFT = fft(cleanFrame);
        noisyFFT = fft(noisyFrame);

        % PSD estimation
        Sxx(:, i) = Sxx(:, i) + abs(cleanFFT).^2 / nrmics;
        noiseFFT = noisyFFT - cleanFFT;
        Nxx(:, i) = Nxx(:, i) + abs(noiseFFT).^2 / nrmics;
    end
end
%%
% Initialize W
W = complex(zeros(K, L, nrmics));

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
        stackedS(k, l) = sum(conj(W(k, l, :)) .* fftResult(k, l, :))/16;
    end
end


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
