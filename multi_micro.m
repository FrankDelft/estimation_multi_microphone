clear;close all;clc;
load("Data.mat");
fs=16000;
%calculate necessary variables for overlap procdure
window_length=(20*10^-3)*fs;
window=hann(window_length);

overlap=0.5;
N=length(Data(:,1));
D=floor(window_length*(1-overlap));
%K=number of frames
K=floor((N-window_length+D)/D);
%% 
%calculate the DFT of overlapped the clean speech signal
S_clean=zeros(K,window_length);
for i = 1:K
        % Extract segment
        start_frame=1 + (i - 1) * (D);
        end_frame=(i-1) * (D)+window_length;
        frame_data = Clean(start_frame : end_frame,1);
        % Compute dft coefficeints
        fft_segment = fft(frame_data .* window);
        S_clean(i,:)=fft_segment;
end
 


%%
%number of microphones to use
m_num=16;

%split sound signal into separate overlapped frames for each microphone
Y_l_k=zeros(K,window_length);
%array containing all windows for all microphones
Y_m_l_k=zeros(m_num,K,window_length);
var_emp=zeros(m_num,1);
%iterate over all microphones
for m =1:m_num
    %iterate over K frames
    for i = 1:K
        % Extract segment
        start_frame=1 + (i - 1) * (D);
        end_frame=(i-1) * (D)+window_length;
        frame_data = Data(start_frame : end_frame,m);
        % Compute dft coefficeints
        fft_segment = fft(frame_data .* window);
        Y_l_k(i,:)=fft_segment;
    end
Y_m_l_k(m,:,:)=Y_l_k;

%now average over the m microphones
Y_ave_l_k=squeeze(mean(Y_m_l_k(1:m,:,:),1));
%calculate the Variance
var_emp(m)=(1/window_length)*(1/K)*sum(abs(Y_ave_l_k- S_clean).^2,"all");
end


Noise_duration=1;
L_1s=(Noise_duration-0.02)/(0.02*(0.5))+1;
estimated_noise=zeros(m_num,1);

for m = 1:m_num
    mic_m_sum=0;
    for i = 1:L_1s
        start_frame=1 + (i - 1) * (D);
        end_frame=(i-1) * (D)+window_length

        noise_FFT=fft(Data(start_frame:end_frame,m).*window,K);
        mic_m_sum=mic_m_sum+noise_FFT.*(1/L_1s);
    end
    estimated_noise(m)=var(mic_m_sum);
end

%now lets get the fischer information for different numbers of microphones
cum_sum_inverse_estimated_noise=cumsum(1./estimated_noise);

%now lets calculate the CRLB for different numbers of microphones
crlb_m=1./cum_sum_inverse_estimated_noise;

%% 

%reconstruct signal
y_reconstructed=zeros(N,1);
%iterate over K frames
for i = 1:K
    inverse_fourier=ifft(Y_ave_l_k(i,:));
    start_frame=1 + (i - 1) * (D);
    end_frame=(i-1) * (D)+window_length;
    y_reconstructed(start_frame:end_frame)=y_reconstructed(start_frame:end_frame)+inverse_fourier';
end

audiowrite("recponstructed.wav", y_reconstructed, 16000);
audiowrite("m1.wav", Data(:,1), 16000);
%% 

plot(crlb_m)
hold on;
plot(var_emp)
hold off;
legend("CRLB", "Actual")








