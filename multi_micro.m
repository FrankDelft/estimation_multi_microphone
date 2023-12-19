clear;close all;clc;
load("Data.mat");
fs=16000;
%calculate necessary variables for overlap procdure
window_length=(20*10^-3)*fs;
window=hann(window_length);
%calculate energy of the window
U=sum(window.^2).*1/length(window);

overlap=0.5;
N=length(Data(:,1));
D=floor(window_length*(1-overlap));
%K=number of frames
K=floor((N-window_length+D)/D);
%% 
%get an estimate for the noise from first second of 16 microphones, (we could just do periodogram averaging over frames and microphones)



%%
%split sound signal into separate overlapped frames for each microphone
Y_l_k=zeros(K,window_length);
%
Y_m_l_k=zeros(16,K,window_length);

%iterate over all microphones
for m =1:16
    %iterate over K frames
    for i = 1:K
        % Extract segment
        start_frame=1 + (i - 1) * (D);
        end_frame=(i-1) * (D)+window_length;
        frame_data = Data(start_frame : end_frame,m);
        % Compute periodogram and store phase information
        fft_segment = fft(frame_data .* window);
        Y_l_k(i,:)=fft_segment;
    end
Y_m_l_k(m,:,:)=Y_l_k;
end





%% 

% Periodogram averaging Welch
function PSD = welch(audioData,L,Overlap,window)
    % percentage overlap
    N=length(audioData);
    D=L*(1-Overlap);
    K=floor((N-L+D)/D);
    segmented_data = zeros(L, K);
    U=sum(window.^2).*1/length(window);
    
    P_i_hat_w = zeros(L, K);
    for i = 1:K
        segmented_data(:, i) = audioData(1 + (i - 1) * (D) : (i-1) * (D)+L);
        P_i_hat_w(:, i) = (1/U).*(1 / L) .* abs(fft(segmented_data(:, i).*window)).^2;
    end
    P_ave_w = zeros(1, L);
    for i = 1:K
        P_ave_w(:) = P_ave_w(:) + P_i_hat_w(:, i);
    end
    PSD = P_ave_w .* (1 / K);
end