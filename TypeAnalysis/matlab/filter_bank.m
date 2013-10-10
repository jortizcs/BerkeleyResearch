%% filter bank implementation in matlab. 
% You can edit this code for any sampling and cutoff frequencies by simply changing the corresponding values. Read the comments and just make the appropriate edits
% Author: Avinash Parnandi,  http://robotics.usc.edu/~parnandi/

%here data is the signal that will be passed through the filter bank
raw_data = load('/home/romain/Unison/Exp/GUTP/dataBerkeley/CoryHallSorted2011/CoryHall_3pw_elt-A_AF_Elevator_.dat');
data = test(:,2);

f=1/0.3333333; %% sampling frequency in MINUTES (we have one point every 20 seconds = 0.3333333 minutes)

order = 3; % order of the filter is 3


fnorm1 = [1/20]/(f/2); % for bandpass, here 1 and 3 are the lower and upper cutoff respectively
% [b,a] = butter(n,Wn) % n = order of the filter, wn is the normalized
% cutoff freq
[b1,a1] = butter(order,fnorm1,'high'); 
data_1 = filtfilt(b1,a1,data); % band pass filter

fnorm2 = [1/360 1/20]/(f/2); % for bandpass
[b2,a2] = butter(order,fnorm2);
data_2 = filtfilt(b2,a2,data); % band pass filter

fnorm3 = [1/360]/(f/2); % for bandpass
[b3,a3] = butter(order,fnorm3,'low');
data_3 = filtfilt(b3,a3,data); % band pass filter


data_4 = data - (data_1+data_2+data_3); % residual

figure
 subplot(5,1,1), plot(data), title('Original data')
 subplot(5,1,2), plot(data_1), title('Highpass filter cutoff at 1/20minutes')
 subplot(5,1,3), plot(data_2), title('Bandpass filter (1/20minutes - 1/360minutes)')
 subplot(5,1,4), plot(data_3), title('Lowpass filter cutoff at 1/360minutes')
 subplot(5,1,5), plot(data_4), title('residual = original data - filters results')
 

