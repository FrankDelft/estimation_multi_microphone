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
%number of microphones to use
m_num=2;

%split sound signal into separate overlapped frames for each microphone
Y_l_k=zeros(K,window_length);
%array containing all windows for all microphones
Y_m_l_k=zeros(m_num,K,window_length);

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
end

%now average over the 16 microphones
Y_ave_l_k=squeeze(mean(Y_m_l_k,1));

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

%the time averaged gives pretty much the same thing as the  
y_time_ave=squeeze(mean(Data,2));
% Subplot with two plots horizontally
subplot(1, 2, 1); % 1 row, 2 columns, first plot
plot(y_reconstructed);
xlabel('Samples');
ylabel('Magnitude');
title('Averaged Y');

subplot(1, 2, 2); % 1 row, 2 columns, second plot
plot(Data(:,1));
xlabel('Samples');
ylabel('Magnitude');
title('Original y');
linkaxes([subplot(1, 2, 1), subplot(1, 2, 2)], 'y');




%% 
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
%calculate the Variance
%1st microphone variance 
var_emp_m1_S=(1/window_length)*(1/K)*sum(abs((squeeze(Y_m_l_k(1,:,:)))- S_clean).^2,"all");

var_emp=(1/window_length)*(1/K)*sum(abs(Y_ave_l_k- S_clean).^2,"all");




