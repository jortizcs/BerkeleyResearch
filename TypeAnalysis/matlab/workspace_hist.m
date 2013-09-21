addpath('/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/matlab')


sz = 20000;
base = '~/StreamFS/StreamFS.classifier/data/KETI/';
cd (base);
folders = dir(base);
DATA=[];
bins = 10;

% for ground truth
PIR=3;
CO2=1;
TEMP=4;
HUM=2;
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
        type='';
        if ~isempty(strfind(fname, 'co2'))
            GT=[GT; CO2];
            type='CO2';
        elseif ~isempty(strfind(fname, 'humidity'))
            GT=[GT; HUM];
            type='HUM';
        elseif ~isempty(strfind(fname, 'light'))
            GT=[GT; PIR];
            type='PIR';
        elseif ~isempty(strfind(fname, 'temperature'))
            GT=[GT; TEMP];
            type='TEMP';
        end
        fprintf('\tGround Truth  [%s]\n',type);
        
        fprintf('\tOpening: %s, %d entries\n', strcat(thisdir, fname), sz);
        fid=fopen(strcat(thisdir, fname));
        vec = fscanf(fid,'%f',[sz,1]);
        fclose(fid);
        
        [N,X]=hist(vec,bins);
        X= sort(X);
        DATA=[DATA; X(bins) X(bins-1) X(bins-2)];
        fprintf('\tsize(DATA)=[%d %d]\n\n', size(DATA,1), size(DATA,2));
    end
    %fprintf('\n');
end

fprintf('Done\n\n\n');



%% Eval success
mid = floor(size(DATA,1)/2);
len = length(DATA);
mdl = ClassificationKNN.fit(DATA(1:mid,:),GT(1:mid));

cnt=0;
for i=mid+1:len
   if predict(mdl,DATA(i,:))==GT(i)
       cnt=cnt+1;
   end
end

cnt/(len-mid)


%% PLOT
figure;
for i=1:size(GT,1)
    if GT(i,1)==PIR
        plot3(DATA(i,1), DATA(i,2), DATA(i,3),'bo');
    elseif GT(i,1)==CO2
        plot3(DATA(i,1), DATA(i,2), DATA(i,3),'g*');
    elseif GT(i,1)==TEMP
        plot3(DATA(i,1), DATA(i,2), DATA(i,3),'rx');
    elseif GT(i,1)==HUM
        plot3(DATA(i,1), DATA(i,2), DATA(i,3),'cx');
    end
    hold on;
end
hold off;
























