clear all;close all;clc
%load data
load('Data.mat'); % Assuming Data.mat contains 'Data' and 'Clean'
% Parameters
fs = 16000; % Sampling frequency 
nrmics=16;
varempEst_wiener = zeros(nrmics, 1);
varempEst_ave = zeros(nrmics, 1);

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



for nrmicsManual = 1:nrmics
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
  
    %variance over L frames
    [var_m, estimated_noise_psd]=noiseEstimation(nrmicsManual,K,frameLength,Data,window,shiftSize);
    % Estimate PSDs for signal and noise
    Sxx = zeros(K, L);
    Nxx = zeros(K, L);
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

    % Initialize W
    W = complex(zeros(K, L, nrmicsManual));
    % Calculate Wiener filter weights
    for k = 1:K
        for l = 1:L
            % Wiener filter
            W(k, l, :) = Sxx(k, l) ./ (Sxx(k, l) + Nxx(k, l));
        end
    end

    % Perform the multiplication and sum
    for k = 1:K
        for l = 1:L
            % Element-wise multiplication and sum across the M dimension
            stackedS(k, l) = sum(conj(W(k, l, :)) .* fftResult(k, l, :))/nrmicsManual;
        end
    end
    
    %calculate empirical Variance
    varempEst_wiener(nrmicsManual)=empVarCalc(K,L,shiftSize,Clean, stackedS,window);

    %now lets do a averaging of the microphones to see the performance
    S_mic_ave=mean(fftResult,3);
    varempEst_ave(nrmicsManual)=empVarCalc(K,L,shiftSize,Clean,S_mic_ave ,window);
    
end
%% 

%now lets get the fischer information for different numbers of microphones
cum_sum_inverse_estimated_noise=cumsum(1./var_m);

%now lets calculate the CRLB for different numbers of microphones
crlb_m=1./cum_sum_inverse_estimated_noise;

plot(crlb_m)
hold on;
plot(varempEst_wiener)
hold off;
%legend("crlb","actual")
   
reconstructedSignal_wiener=reconstructTimeDomainSignal(stepSize,frameSize,numFrames,stackedS);
reconstructedSignal_ave=reconstructTimeDomainSignal(stepSize,frameSize,numFrames,S_mic_ave);
figure;
plot(Clean, 'b', 'DisplayName', 'Original');
hold on;
plot(reconstructedSignal_wiener, 'r', 'DisplayName', 'Modeled');

% Calculate the Mean Squared Error (MSE)
mseValue_wiener = mean((Clean((1:end-55)) - reconstructedSignal_wiener).^2); % some part of the reconstructedSignal (last 55 samples) was lost due to calculations so we had to exclude the last 55 samples of the clean signal 
mseValue_ave = mean((Clean((1:end-55)) - reconstructedSignal_ave).^2); 
mseValue_m1=mean((Clean((1:end-55)) - Data(1:end-55,1)).^2);%the first microphones MSE value

xlabel('Time');
ylabel('Value');
title('Original vs Modeled Time Series');
legend('show');
grid on;

function [var_freq_ave, var_m_k]= noiseEstimation(nrmics,K,frameLength,Data,window,shiftSize)
    %now to make an estimate of the noise based on the 1st second of noise
    Noise_duration=1;
    L_1s=(Noise_duration-frameLength)/(frameLength*(0.5))+1;
    %variance over L frames for fft (same as psd)
    estimated_noise_psd=zeros(nrmics,K);

    for m = 1:nrmics
        mic_m_sum=zeros(K,1);
        for l = 1:L_1s
            frameStart = (l-1) * shiftSize + 1;
            frameEnd = frameStart + K - 1;
            
            noisyFrame=Data(frameStart:frameEnd,m).*window;
            noise_FFT=fft(noisyFrame,K);
            mic_m_sum=mic_m_sum+abs(noise_FFT).^2;
        end
        estimated_noise_psd(m,:)=(mic_m_sum').*1/L_1s;
    end
    %average over the frequency bins and 
    var_freq_ave=mean(estimated_noise_psd,2);
    var_m_k=estimated_noise_psd;
end

%reconstruct time domain signal
function reconstructedSignal=reconstructTimeDomainSignal(stepSize,frameSize,numFrames,fftSignal)
    reconstructedSignal = zeros((numFrames-1) * stepSize + frameSize, 1);
    
    % Overlap-add method
    for i = 1:numFrames
        frameStart = (i-1) * stepSize + 1;
        frameEnd = frameStart + frameSize - 1;
        % Inverse FFT
        ifftFrame = ifft(fftSignal(:, i), 'symmetric');
        % Overlap-add
        reconstructedSignal(frameStart:frameEnd) = reconstructedSignal(frameStart:frameEnd) + ifftFrame;
    end
    audiowrite("reconstructed.wav", reconstructedSignal, 16000);
end


function varEmpEst = empVarCalc(K,L,shiftSize,Clean, signalEst,window )
    varianceSum = 0;
    % Process each frame
    for l = 1:L
        frameStart = (l-1) * shiftSize + 1;
        frameEnd = frameStart + K - 1;
        % Frame from clean signal
        cleanFrame = Clean(frameStart:frameEnd);
        % FFT of frames
        S = fft(cleanFrame.*window, K);
        % Variance for this frame
        varianceSum = varianceSum + sum(abs(signalEst(:,l) - S).^2);
    end
        % Average variance for this number of microphones
        varEmpEst = varianceSum / (K * L);
end
