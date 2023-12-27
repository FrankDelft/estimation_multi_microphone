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
    
w = complex(rand(nrmics, 1), rand(nrmics, 1)); % Complex noise reduction weights
    
for i = 1:numFrames
    % Initialize a temporary variable to sum the weighted FFT results for this frame
    combinedFrame = zeros(frameSize, 1);
    for m = 1:nrmics
        % Multiply FFT result by the weight for this microphone
        weightedFrame = fftResult(:, i, m) * w(m);

        % Sum the weighted results across microphones
        combinedFrame = combinedFrame + weightedFrame;
    end
    % Store the combined result for this frame
    stackedS(:, i) = combinedFrame;
end

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

plot(reconstructedSignal)
