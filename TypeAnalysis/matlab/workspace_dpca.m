%% Code to process the data using EMD + DPCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
cd '~/StreamFS/StreamFS.classifier/';

%%
files = dir ('data/KETI/413');
dims = size(files,1);
% figure;

for i=1:dims
    fn = files(i,1);
    
    if ~isempty(strfind(fn.name, 'temperature')) || ~isempty(strfind(fn.name, 'co2'))
        fprintf('processing: %s\n', fn.name);
        f = importdata(strcat('data/KETI/413/', fn.name));
%         plot(f);
%         title(fn.name);
    end
end
%clear f;
clear dims;
clear files;

%% 
% subf = f(1:10000);
% plot(subf);
cd '~/StreamFS/StreamFS.classifier/';

f1 = importdata('data/KETI/413/temperature.csv');
f1sub = f1(1:20000,1);
clear f1;
f2 = importdata('data/KETI/413/co2.csv');
f2sub = f2(1:20000,1);
clear f2;
f3 = importdata('data/KETI/413/light.csv');
f3sub = f3(1:20000,1);
clear f3;
f4 = importdata('data/KETI/424/co2.csv');
f4sub = f4(1:2000,1);
clear f4;


% EMD


% IMF_f1 = emd(f1sub);%,'MAXMODES',3);
% IMF_f2 = emd(f2sub);%,'MAXMODES',3);
% IMF_f3 = emd(f3sub);%,'MAXMODES',3);

%%  re-aggregate IMFs
addpath('/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/matlab')

[IMF_f1, REAGG_IMF_f1] = StripAgg(f1sub,1/15);

% figure;
% for i=1:size(REAGG_IMF_f1,1)
%     subplot(size(REAGG_IMF_f1,1),1,i);
%     plot(REAGG_IMF_f1(i,:),'b');
% end
% 
% figure;
% for i=1:size(IMF_f1,1)
%     subplot(size(IMF_f1,1),1,i);
%     plot(IMF_f1(i,:),'color','b');
% end





[IMF_f2, REAGG_IMF_f2] = StripAgg(f2sub,1/15);

% figure;
% for i=1:size(REAGG_IMF_f2,1)
%     subplot(size(REAGG_IMF_f2,1),1,i);
%     plot(REAGG_IMF_f2(i,:),'b');
% end

% figure;
% for i=1:size(IMF_f2,1)
%     subplot(size(IMF_f2,1),1,i);
%     plot(IMF_f2(i,:),'color','b');
% end






[IMF_f3, REAGG_IMF_f3] = StripAgg(f3sub,1/15);

% figure;
% for i=1:size(REAGG_IMF_f3,1)
%     subplot(size(REAGG_IMF_f3,1),1,i);
%     plot(REAGG_IMF_f3(i,:),'b');
% end

% figure;
% for i=1:size(IMF_f3,1)
%     subplot(size(IMF_f3,1),1,i);
%     plot(IMF_f3(i,:),'color','b');
% end


[IMF_f4, REAGG_IMF_f4] = StripAgg(f4sub,1/15);

% figure;
% for i=1:size(REAGG_IMF_f4,1)
%     subplot(size(REAGG_IMF_f4,1),1,i);
%     plot(REAGG_IMF_f4(i,:),'b');
% end

%%
clear all;
f4_temp =importdata('~/StreamFS/StreamFS.classifier/data/KETI/424/temperature.csv');
f4sub = f4_temp(1:20000,1);
clear f4_temp;
[IMF, REAGG] = StripAgg(f4sub,1/15);

figure;
for i=1:size(REAGG,1)
    subplot(size(REAGG,1),1,i);
    plot(REAGG(i,:),'b');
end


%%
LF_IMF = 4;
MF_IMF = 3;
HF_IMF = 2;
RAW_IMF = 1;
RES = 5;

IMF_band = RES; 

if IMF_band==LF_IMF
    band_name = ' LF IMF';
elseif IMF_band==MF_IMF
    band_name = ' MF IMF';
elseif IMF_band == HF_IMF;
    band_name = ' HF IMF';
elseif IMF_band == RES
    band_name = ' Residual';
elseif IMF_band == RAW_IMF
    band_name = ' Raw Data';
end

figure;
subplot(4,1,1);
plot(REAGG_IMF_f1(IMF_band,:),'color','b');
title(strcat('temperature ',band_name));

subplot(4,1,2);
plot(REAGG_IMF_f2(IMF_band,:),'color','r');
title(strcat('co2 ', band_name));

subplot(4,1,3);
plot(REAGG_IMF_f3(IMF_band,:),'color','g');
title(strcat('light ', band_name));

subplot(4,1,4);
plot(REAGG_IMF_f4(IMF_band,:),'color','y');
title(strcat('co2_2 ', band_name));

%%

% separate it out
w = 100;

F1 =[];
for i=0: (size(REAGG_IMF_f1',1)/w)-1
    start_w = i*w+1;
    end_w = start_w+w-1;
    F1=[F1; REAGG_IMF_f1(IMF_band,start_w:end_w)];
end

F2 =[];
for i=0: (size(REAGG_IMF_f2',1)/w)-1
    start_w = i*w+1;
    end_w = start_w+w-1;
    F2=[F2; REAGG_IMF_f2(IMF_band,start_w:end_w)];
end

F3 =[];
for i=0: (size(REAGG_IMF_f3',1)/w)-1
    start_w = i*w+1;
    end_w = start_w+w-1;
    F3=[F3; REAGG_IMF_f3(IMF_band,start_w:end_w)];
end

F4 =[];
for i=0: (size(REAGG_IMF_f4',1)/w)-1
    start_w = i*w+1;
    end_w = start_w+w-1;
    F4=[F4; REAGG_IMF_f4(IMF_band,start_w:end_w)];
end


% PCA
DATA = [F1;F2;F3;F4];


start_ = 1;
end_= size(F1,1);
[wcoeff,score,latent,tsquared,explained] = pca(DATA, 'VariableWeights','variance');
figure; plot3(score(start_:end_,1), score(start_:end_,2), score(start_:end_,3), 'bo');

start_ = end_+1;
end_= start_+size(F2,1)-1;
hold on; plot3(score(start_:end_,1), score(start_:end_,2), score(start_:end_,3), 'rx');

start_ = end_+1;
end_= start_+size(F3,1)-1;
hold on; plot3(score(start_:end_,1), score(start_:end_,2), score(start_:end_,3), 'g*');

start_ = end_+1;
end_= start_+size(F4,1)-1;
hold on; plot3(score(start_:end_,1), score(start_:end_,2), score(start_:end_,3), 'y*')
grid on

%% Ground truth

TEMP=3;
CO2=2;
LIGHT=1;

start_ = 1;
end_= size(F1,1);
GT=[ones(end_-start_+1,1)*TEMP];

start_ = end_+1;
end_= start_+size(F2,1)-1;
GT=[GT; ones(end_-start_+1,1)*CO2];

start_ = end_+1;
end_= start_+size(F3,1)-1;
GT=[GT;ones(end_-start_+1,1)*LIGHT];

start_ = end_+1;
end_= start_+size(F4,1)-1;
GT=[GT;ones(end_-start_+1,1)*CO2];

%% kmeans
IDX= kmeans(score(:,1),3);

%%
same_cnt=0;
figure;
for i=1:length(IDX)
    if IDX(i)==1
        scatter(score(i),0, 'g*');
    elseif IDX(i)==2
        scatter(score(i),0, 'r*');
    elseif IDX(i)==3
        scatter(score(i),0, 'b*');
    end
    hold on;
    
    if IDX(i)==GT(i)
        same_cnt=same_cnt+1;
    end
end
hold off;

hit_rate=same_cnt/length(IDX)





%% Full analysis
addpath('/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/matlab')

% IMF bands
LF_IMF = 4;
MF_IMF = 3;
HF_IMF = 2;
RAW_IMF = 1;
RES = 5;
IMF_band = RES; 


sz = 10000;
base = '~/StreamFS/StreamFS.classifier/data/KETI/';
cd (base);
folders = dir(base);
srate = 1/15; %sample rate
w = 100;    %window for dpca
DATA=[];

% for ground truth
PIR=3;
CO2=1;
TEMP=2;
HUM=4;
GT=[];

start_=0;
end_=0;

for i=4:size(folders,1)-1
    fname= folders(i).name;
    files = dir(strcat(base,fname));
    thisdir = strcat(strcat(base, fname),'/');
    fprintf('Processing %s.\n',thisdir);

    for j=3:size(files,1)
        fname = files(j).name;
        
        % record the type for ground truth
        start_ = end_+1;
        end_ = start_-1+(sz/w);
        type='';
        if ~isempty(strfind(fname, 'co2'))
            GT=[GT; CO2, start_, end_];
            type='CO2';
        elseif ~isempty(strfind(fname, 'humidity'))
            GT=[GT; HUM, start_, end_];
            type='HUM';
        elseif ~isempty(strfind(fname, 'light'))
            GT=[GT; PIR, start_, end_];
            type='PIR';
        elseif ~isempty(strfind(fname, 'temperature'))
            GT=[GT; TEMP, start_, end_];
            type='TEMP';
        end
        fprintf('\tGround Truth  [%s %d %d]\n',type,start_, end_);
        
        fprintf('\tOpening: %s, %d entries\n', strcat(thisdir, fname), sz);
        fid=fopen(strcat(thisdir, fname));
        vec = fscanf(fid,'%f',[sz,1]);
        fclose(fid);
        
        % EMD dis-aggregation
        fprintf('\tEMD Running\n');
        [IMFS, REAGG] = StripAgg(vec,srate);
        
        
        % separate it out into segments of size w
        fprintf('\tPartitioning, window=%d\n',w);
        SIG =[];
        for k=0: (size(REAGG',1)/w)-1
            start_w = k*w+1;
            end_w = start_w+w-1;
            SIG=[SIG; REAGG(IMF_band,start_w:end_w)];
        end
        DATA=[DATA; SIG];
        fprintf('\tsize(DATA)=[%d %d]\n\n', size(DATA,1), size(DATA,2));
    end
    %fprintf('\n');
end

fprintf('Done\n\n\n');

%%
% PCA
[wcoeff,score,latent,tsquared,explained] = pca(DATA, 'VariableWeights','variance');


%%
score2 = zeros(size(score,1), 3);
for i=1:length(score2)
    score2(i,1)=score(i,1);
    score2(i,2)=score(i,2);
    score2(i,3)=score(i,1)^2 + score(i,2)^2;
end
%IDX= kmeans(score(:,1),4);
%IDX= kmeans(score2,4);

figure;
for i=1:size(GT,1)
    start_ = GT(i,2);
    end_ = GT(i,3);
    if GT(i,1)==PIR
        plot3(score2(start_:end_,1), score2(start_:end_,2), score2(start_:end_,3), 'bo');
    elseif GT(i,1)==CO2
        plot3(score2(start_:end_,1), score2(start_:end_,2), score2(start_:end_,3), 'g*');
    elseif GT(i,1)==TEMP
        plot3(score2(start_:end_,1), score2(start_:end_,2), score2(start_:end_,3), 'rx');
    elseif GT(i,1)==HUM
        plot3(score2(start_:end_,1), score2(start_:end_,2), score2(start_:end_,3), 'cx');
    end
    hold on;
end
hold off;

grid on;



%%
% Plot
for i=1:size(GT,1)
    start_ = GT(i,2);
    end_ = GT(i,3);
    if GT(i,1)==PIR
        plot3(score(start_:end_,1), score(start_:end_,2), score(start_:end_,3), 'bo');
    elseif GT(i,1)==CO2
        plot3(score(start_:end_,1), score(start_:end_,2), score(start_:end_,3), 'g*');
    elseif GT(i,1)==TEMP
        plot3(score(start_:end_,1), score(start_:end_,2), score(start_:end_,3), 'rx');
    elseif GT(i,1)==HUM
        plot3(score(start_:end_,1), score(start_:end_,2), score(start_:end_,3), 'cx');
    end
    hold on;
end
hold off;




%%  EVALUATION
NGT=zeros(GT(size(GT,1),3),1);
for i=1:size(GT,1)
    for j=GT(i,2):GT(i,3)
        NGT(j,1)=GT(i,1);
    end
end

% Eval success
mid = floor(size(DATA,1)/2);
len = length(DATA);
mdl = ClassificationKNN.fit(DATA(1:mid,1:3),NGT(1:mid));

cnt=0;
for i=mid+1:len
   if predict(mdl,DATA(i,1:3))==NGT(i)
       cnt=cnt+1;
   end
end

cnt/(len-mid)

% COMP = [NGT, IDX];
% 
% length(IDX(find(IDX(:,1)~=NGT(:,1))))/length(IDX)

% RES=58.8%
% RAW = 70.84% w/o energy & 76.5 w/energy%

%%
total=[];
for j=1:4 %classes (labled)
    r = find(GT(:,1)==j);
    SUBGT = [GT(r,2), GT(r,3)];
    SUBIDX = [];
    for b=1:size(SUBGT,1)
        SUBIDX=[SUBIDX; IDX(SUBGT(b,1):SUBGT(b,1))];
    end
    
    per1 = length(find(SUBIDX(:,1)==1))/length(SUBIDX);
    per2 = length(find(SUBIDX(:,1)==2))/length(SUBIDX);
    per3 = length(find(SUBIDX(:,1)==3))/length(SUBIDX);
    per4 = length(find(SUBIDX(:,1)==4))/length(SUBIDX);
    if j==PIR
        t='PIR';
    elseif j==CO2
        t='CO2';
    elseif j==TEMP
        t='TEMP';
    elseif j==HUM
        t='HUM';
    end
    fprintf('[%s %f %f %f %f]\n', t, per1, per2, per3,per4);
end

%%
G1=length(find(GT(:,1)==1));
G2=length(find(GT(:,1)==2));
G3=length(find(GT(:,1)==3));
G4=length(find(GT(:,1)==4));
fprintf('[GT %d %d %d %d]\n',G1,G2,G3,G4); 


I1=length(find(IDX(:,1)==1));
I2=length(find(IDX(:,1)==2));
I3=length(find(IDX(:,1)==3));
I4=length(find(IDX(:,1)==4));
fprintf('[IDX %d %d %d %d]\n',G1,G2,G3,G4);














