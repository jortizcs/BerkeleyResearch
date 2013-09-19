function [ IMF, REAGG ] = StripAgg( s, sampling_rate )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    IMF = emd(s);
    
    HF_MAX = 1/(20*60);         % 20 minutes
    MF_MAX = 1/(60*60*6);       % 6 hours
    LF_MAX = 1/(60*60*24*6);    % 6 days

    %REAGG = [];
    LAGG=zeros(1,size(IMF,2));
    MAGG=zeros(1,size(IMF,2));
    HAGG=zeros(1,size(IMF,2));
    RAGG=zeros(1,size(IMF,2));

    if isempty(sampling_rate)
        sampling_rate = 1;
    end

    for i=1:size(IMF,1)
        avgT = ZCR(IMF(i,:))*sampling_rate;
        if avgT>HF_MAX
            HAGG = HAGG + IMF(i,:);
        elseif avgT<=HF_MAX && avgT>MF_MAX
            MAGG = MAGG + IMF(i,:);
        elseif avgT<=MF_MAX && avgT>LF_MAX
            LAGG = LAGG + IMF(i,:);
        else
            RAGG = RAGG + IMF(i,:);
        end
    end
    REAGG = [s';HAGG;MAGG;LAGG;RAGG];

end

