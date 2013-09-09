clear
clc
%[audioData, fs, bps] = wavread('C:\Users\dh5gm\Desktop\Data\sound clip\STG_LungS_Norm_Tracheal.wav');
[audioData, fs, bps] = wavread('C:\Users\dh5gm\Desktop\Data\sound clip\STG_LungS_Asthma_Wz.wav');
%[audioData, fs, bps] = wavread('C:\Users\dh5gm\Desktop\sound clip\STG_LungS_Pneumonia_Cr.wav');
%[audioData, fs, bps] = wavread('C:\Users\dh5gm\Desktop\sound clip\STG_LungS_CHF_Crackles.wav');
%[audioData, fs, bps] = wavread('C:\Users\dh5gm\Desktop\sound clip\normal_tracheal.wav');
%[audioData, fs, bps] = wavread('C:\Users\dh5gm\Desktop\sound clip\wheezing_a.wav');
%[audioData, fs, bps] = wavread('C:\Users\dh5gm\Desktop\sound clip\right2.wav');
audioData = audioData(:,1);
length = size(audioData,1);%number of samplings
t = length/fs*linspace(0,1,length);%time as x-axis

%filter
w = window(@blackman, floor(fs/12));
B = fir1(floor(fs/12)-1, [100,2000]/(fs/2), w); %window-based bandpassFIR
%B = fir1(floor(fs/20)-1, 2/(fs/2), 'low', w); %window-based FIR
A = 1;
newData = filter(B, A, audioData);
x = newData;
% wn = [50, 1000]/(fs/2);
% [B, A] = butter(1, wn, 'stop');
% newData = filter(B, A, AudioData);
% x = newData;  
% newData = abs(newData); %for computing the cpm

% %downsample by 100
% dsData = downsample(newData, 50);
% dsT = downsample(t, 50);

%average, downsample by 100
% j = 1;
% aveData = zeros(length/100,1);
% for k= 1:size(aveData,1)*100
%     if mod(k,100)~=0
%         aveData(j) = aveData(j) + newData(k);
%     else
%         aveData(j) = aveData(j) + newData(k);
%         aveData(j) = aveData(j)/100;
%         j = j+1;
%     end
% end

% %de-noise
% length = size(dsData,1);
% dnData = zeros(length,1);
% tmp = 0;
% for k= 1:length
%     if mod(k,80)~=0
%         tmp = tmp + dsData(k);
%     else
%         tmp = tmp + dsData(k);
%         tmp = tmp/80;
%         for j = k-79:k
%             if j == 1
%                 dnData(j) = 0;
%             elseif dsData(j)>2*tmp
%                 dnData(j) = dnData(j-1);
%             else
%                 dnData(j) = dsData(j); 
%             end
%         end
%         tmp = 0;
%     end
% end

% %downsample by 4
% dssData = downsample(dnData, 4);
% dsTT = downsample(t,200);

%average, downsample again by 10
% length = size(dnData,1);
% j = 1;
% dnspData = zeros(length/10,1);
% for k= 1:size(dnspData,1)*10
%     if mod(k,10)~=0
%         dnspData(j) = dnspData(j) + dnData(k);
%     else
%         dnspData(j) = dnspData(j) + dnData(k);
%         dnspData(j) = dnspData(j)/10; 
%         j = j+1;
%     end
% end

%LPF
% w = window(@blackman, fs/200*2);
% B = fir1(fs/200*2-1, 0.6/(fs/2/200), 'low', w); %window-based FIR
% A = 1;
% lpfData = filter(B, A, dssData);
% wn = 0.6/(fs/2/200);
% [B, A] = butter(3, wn, 'low');
% lpfData = filter(B, A, dssData);

%derivative
% dt = diff(lpfData);
% dt = dt/(dsTT(2)-dsTT(1));

%spectrogram
figure
% spectrogram(x,256,200,256,fs)
[S,F,T,P] = spectrogram(x,256,200,256,fs);
surf(T,F,10*log10(abs(P)),'EdgeColor','none');axis tight; 
view(0,90)
xlabel('Time (Seconds)')
ylabel('Frequency (Hz)')

%PSD
%x = newData;  
%x = AudioData;
% h = spectrum.welch;                  % Create a Welch spectral estimator. 
tt = 0.1;
peak = zeros(10,3);
ctr = 1;
figure %for distribution
hold on %for distribution
while tt<max(t)
    str = sprintf('Start T = %.2fs', tt);
    % figure %for psd
    % [psd1, F] = pyulear(x(tt*fs:(tt+0.3)*fs),4,fs);
    % [psd2, F] = pyulear(x(tt*fs:(tt+0.3)*fs),25,fs);
    % [psd3, F] = pyulear(x(tt*fs:(tt+0.3)*fs),50,fs);
    span = 6;
    for k = 1:span
        [psd, F] = pyulear(x(tt*fs:(tt+0.3)*fs),100,fs);
%         subplot(span,1,k) %for psd
%         plot(F/pi*fs/2, psd)  %for psd
%         axis([0 1000 0 1.25*max(psd)])    %for psd
        tt = tt + 0.15;

        %peak & BW
        pos = find(psd == max(psd));%index of peak in psd matrix
        peak(ctr,1) = F(pos);%freq of peak
        low=0;
        high=0;
        tmp = pos;
        while psd(tmp)>0.05*psd(pos)
            tmp = tmp-1;
        end
        low = F(tmp);
        tmp = pos;
        while psd(tmp)>0.05*psd(pos)
            tmp = tmp+1;
        end
        high = F(tmp);
        peak(ctr,2) = high-low;
        ctr = ctr+1;
    end
    plot(peak(:,1), peak(:,2),'r.','MarkerSize',14); %for distribution
    xlabel('Frequency (Hz)')
end
peak = peak/pi*fs/2;
%psd(h,x,'Fs',fs);                    % Calculate and plot the one-sided psd.
%hpsd = psd(h,x,'ConfLevel',.98);     % psd with confidence level
%figure,plot(hpsd)

%cpm
% figure
% subplot(4,1,1)
% plot(t, audioData)
% grid on
% xlabel('Original Data')
% subplot(4,1,2)
% plot(t, newData)
% plot(dsT(1:size(dsData,1)), dsData)
% grid on
% xlabel('Bandpass Filtered and Downsampled')
% subplot(4,1,3)
% plot(dsTT(1:size(dssData,1)), dssData)
% grid on
% xlabel('De-noised and further Downsampled')
% subplot(4,1,4)
% hold on
% plot(dsTT(1:size(lpfData,1)), lpfData)
% plot(dsTT(2:end), dt, 'r:','LineWidth',1.5)
% grid on
% xlabel('LPFed and Derivative')
% lh = legend('signal','derivative');
% set(lh,'Location','SouthWest','Orientation','vertical')
% %plot(t, audioData)