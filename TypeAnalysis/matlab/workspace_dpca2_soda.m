%% Full analysis (slow)
addpath('/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/matlab')

% IMF bands
LF_IMF = 4;
MF_IMF = 3;
HF_IMF = 2;
RAW_IMF = 1;
RES = 5;
IMF_band = RES; 


sz = 5000;
base = '/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/data/';
cd (base);
folders = dir(base);
srate = 1/15; %sample rate
w = 100;    %window for dpca
DATA=[];

start_=0;
end_=0;

%type counts
art=0;
oat=0;
rvav=0;
vav=0;

figure;
max=0;
min=2*10^10;
for i=4:size(folders,1)-1
    fname= folders(i).name;
    files = dir(strcat(base,fname));
    thisdir = strcat(strcat(base, fname));
    %fprintf('Processing %s\n',thisdir);
    
    type = '';
    if ~isempty(strfind(fname, 'ART')) && art==0
        type= 'ART';
        art = 1;
    elseif ~isempty(strfind(fname, 'OAT')) && oat==0
        type = 'OAT';
        oat = 1;
    elseif ~isempty(strfind(fname, 'RVAV')) && rvav==0
        type = 'RVAV';
        rvav=1;
    elseif ~isempty(strfind(fname, 'VAV')) && vav==0
        type = 'VAV';
        vav=1;
    end
    
    if strcmp(type,'')==0
        files = dir(strcat(strcat(base,fname),'/'));
        for j=1:size(files,1)
            if ~isempty(strfind(files(j,1).name,'03M.DAT'))
                fname = strcat(strcat(strcat(base,fname),'/'),files(j,1).name);
                fprintf('\tfname=%s\n',files(j,1).name);
                
                data = zeros(sz,2);
                pts=0;
                fid = fopen(fname);
                line = 0;
                while pts<5000 & line ~= -1
                    %fprintf('pts_read=%d\n',pts);
                    line = fgets(fid);
                    if line ~=-1
                        t=strsplit(line,',\t');
                        ts=t(1,1);
                        v=t(1,6);
                        pts = pts+1;
                        if isempty(strfind(v{1,1},'nan'))
                            data(pts,1)=str2num(ts{1,1});
                            data(pts,2)=str2double(v{1,1});
                            
                            if data(pts,1)>max
                                max=data(pts,1);
                            end
                            if data(pts,1)<min
                                min = data(pts,1);
                            end
                            fprintf('%s range(pts=%d)=[%d %d]\n',files(j,1).name,pts,min,max);
                            
                        end
                    else
                        data = data(find(data(:,1)~=0),:);
                    end
                    
                    DATA = [DATA; data(:,1); data(:,2)];
                end
                fclose(fid);
                
            end
        end
    end
    
    if (art+oat+rvav+vav)==4
        break;
    end
end

%% Full analysis (fast)
addpath('/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/matlab')

% IMF bands
LF_IMF = 4;
MF_IMF = 3;
HF_IMF = 2;
RAW_IMF = 1;
RES = 5;
IMF_band = RES; 


sz = 5000;
base = '/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/data/';
cd (base);
folders = dir(base);
srate = 1/15; %sample rate
w = 100;    %window for dpca
DATA=[];

lengths=zeros(4,1);
types={};

start_=0;
end_=0;

renameAll = '/Users/jortiz/data/rnall2';

%type counts
art=0;
oat=0;
rvav=0;
vav=0;

% figure;
maxts=0;
mints=2*10^10;
for i=4:size(folders,1)-1
    fname= folders(i).name;
    files = dir(strcat(base,fname));
    thisdir = strcat(strcat(base, fname));
    system(strcat(['cp' ' ' renameAll ' ' thisdir]));
    cd(thisdir);
    system('source rnall2');
    
    %fprintf('Processing %s\n',thisdir);
    
    type = '';
    if ~isempty(strfind(fname, 'ART')) && art==0
        type= 'ART';
        art = 1;
    elseif ~isempty(strfind(fname, 'OAT')) && oat==0
        type = 'OAT';
        oat = 1;
    elseif ~isempty(strfind(fname, 'RVAV')) && rvav==0
        type = 'RVAV';
        rvav=1;
    elseif ~isempty(strfind(fname, 'VAV')) && vav==0
        type = 'VAV';
        vav=1;
    end
    
    
    offset = 0;
    len = 5000;
    tlines = 0;
    if strcmp(type,'')==0
        files = dir(strcat(strcat(base,fname),'/'));
        for j=1:size(files,1)
            if ~isempty(strfind(files(j,1).name,'03M.DAT'))
                fname = strcat(strcat(strcat(base,fname),'/'),files(j,1).name);
                [stat, res] = system(strcat(['wc -l' ' ' fname]));
                tok = strtok(res);
                tlines = str2num(tok)
                offset = tlines - offset;
                fprintf('\tfile=%s\n',files(j,1).name);
                
                artcmd = strcat(['~/bin/ack ''(\d+),\t.*,\t(\d+\.\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);
                               
                [stat, res]=system(artcmd);
                artcmd2 = strcat(['~/bin/ack ''(\\d+),\\t.*,\\t(\\d+\\.\\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);
                fprintf(strcat (['artcmd=' artcmd2 '\n']));
                data = importdata('testres.csv');
                system('rm -f testres.csv rnall2');
                DATA = [DATA; data(:,1), data(:,2)];
                
                lengths((art+oat+rvav+vav))=length(data);
                
                if max(data(:,1))>maxts
                    maxts=max(data(:,1));
                end
                if min(data(:,1))<mints
                    mints = min(data(:,1));
                end
                
                types{(art+oat+rvav+vav)}=type;
                
            end
        end
    end
    
    if (art+oat+rvav+vav)==4
        break;
    end
end

%% PLOT THEM
start_ =1
end_ = lengths(1)
for i=1:(art+oat+rvav+vav)
    %plot them
    subplot(4,1,i);
    scatter(DATA(start_:end_,1), DATA(start_:end_,2));
    title(types{i});
    xlim([mints maxts]);
    
    i
    start_ = end_+1
    if i < (art+oat+rvav+vav)
        end_ = start_-1+lengths(i+1)
    end
end

%%
% rnall2 - file.
%
% for file in *.DAT
% do
%     mv "$file" "${file/$/_}"
% done
%%  Classify with PCA (no EMD, DPCA)

start_ =1
end_ = lengths(1)
labels={}
for i=1:(art+oat+rvav+vav)
    %plot them
    subplot(4,1,i);
    
%     scatter(DATA(start_:end_,1), DATA(start_:end_,2));
    %data=[DATA(start_:end_,1), DATA(start_:end_,2)];
    
    for j=1:size(DATA(start_:end_,:),1)
        labels{start_+j-1} = types{i};
    end
    
    start_ = end_+1
    if i < (art+oat+rvav+vav)
        end_ = start_-1+lengths(i+1)
    end
end


%% build model
mdl = ClassificationKNN.fit(DATA,labels');

%% Eval success
addpath('/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/matlab')

% IMF bands
LF_IMF = 4;
MF_IMF = 3;
HF_IMF = 2;
RAW_IMF = 1;
RES = 5;
IMF_band = RES; 


sz = 5000;
base = '/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/data/';
cd (base);
folders = dir(base)
srate = 1/15; %sample rate
w = 100;    %window for dpca
DATA=[];

lengths=zeros(4,1);
class_res={};

start_=0;
end_=0;

renameAll = '/Users/jortiz/data/rnall2';

%type counts
art=0;
oat=0;
rvav=0;
vav=0;

success = 0;
idx =1;

% figure;
maxts=0;
mints=2*10^10;
artcmd='';
for i=4:size(folders,1)-1
    fname= folders(i).name;
    files = dir(strcat(base,fname));
    thisdir = strcat(strcat(base, fname));
    if isdir(thisdir)        
        system(strcat(['cp' ' ' renameAll ' ' thisdir]));
        cd(thisdir);
        system('source rnall2');
   
    
        fprintf('Processing %s\n',thisdir);

        type = '';
        if ~isempty(strfind(fname, 'ART')) && art==0
            type= 'ART';
        elseif ~isempty(strfind(fname, 'OAT')) && oat==0
            type = 'OAT';
        elseif ~isempty(strfind(fname, 'RVAV')) && rvav==0
            type = 'RVAV';
        elseif ~isempty(strfind(fname, 'VAV')) && vav==0
            type = 'VAV';
        end
        
        
        
        if strcmp(type,'')==0
            
            
            
            files = dir(strcat(strcat(base,fname),'/'));
            for j=3:size(files,1)
                offset = 0;
                len = 5000;
                tlines = 0;
                if isempty(strfind(files(j,1).name,'03M.DAT'))
                    
%                     fname = strcat(strcat(strcat(base,fname),'/'),files(j,1).name)
                    fname = strcat(strcat(thisdir,'/'),files(j,1).name);

                    [stat, res] = system(strcat(['wc -l' ' ' fname]));
                    tok = strtok(res);
                    tlines = str2num(tok)
                    offset = tlines-offset
                    fprintf('\tfile=%s\n',files(j,1).name);
                    
                    artcmd = strcat(['~/bin/ack ''(\d+),\t.*,\t(\d+\.\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);
                                   
                    [stat, res]=system(artcmd);
                    artcmd2 = strcat(['~/bin/ack ''(\\d+),\\t.*,\\t(\\d+\\.\\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);
                    fprintf(strcat (['artcmd=' artcmd2 '\n']));
                    data = importdata('testres.csv')
                    system('rm -f testres.csv rnall2');
                    
                    
                    if ~isempty(data)
                        class_res{idx,1} = files(j,1).name;
                        class_res{idx,2} = type;
                        class_res{idx,3} = predict(mdl,data);
    
                        if class_res{idx,2} == class_res{idx,3}
                            success = success +1;
                        end
    
                        idx = idx+1;
                    end
                    
                end
            end
        end
    end
    

end


% success/size(class_res,1)
%% DPCA 

w=100;

start_=1;
end_=0;
SQ_DATA=[];


for i=1:size(lengths,1)
    start_=end_+1;
    end_ = start_-1+lengths(i);
    data = DATA(start_:end_,2);
    nr = floor(size(data,1)/w);
    sq_data = zeros(floor(size(data,1)/w),w);
    
    for j=0:nr-1
        start_w = j*w+1;
        end_w = start_w+w-1;
        sq_data(j+1,:)=sq_data(j+1,:) + data(start_w:end_w,1)';
    end
    SQ_DATA=[SQ_DATA; sq_data];
end

[wcoeff,score,latent,tsquared,explained] = pca(SQ_DATA, 'VariableWeights','variance');

figure;
plot(score(1:7,1),score(1:7,2), 'gx'); %ART
hold on;
plot(score(8:57,1),score(8:57,2), 'ro'); %VAV
hold on;
plot(score(58:107,1),score(58:107,2), 'b*'); %RVAV
hold on;
plot(score(108:157,1),score(108:157,2), 'c+'); %OAT
hold off;

%%
types=cell(157,1);
for i=1:7
    types{i}='ART';
end

for i=8:57
    types{i}='VAV';
end

for i=58:107
    types{i}='RVAV';
end

for i=108:157
    types{i}='OAT';
end

mdl = ClassificationKNN.fit(SQ_DATA,types);
%%

newv=min(data(:,1)):15:max(data(:,1));
u=unique(data(:,1));
n=histc(data(:,1),u);
data([find(n>1)],:)=[]
i1 = interp1(data(:,1),data(:,2),newv', 'spline');

%%
idx =0;
type='OAT';
dicts={'ART',5000;'OAT',5000};
idx=0;
for id=1:size(dicts,1)
    if strcmpi(dicts{id,1},type)==1
        idx=id;
    end
end
idx



% keep track of number learning examples
thisid =0;
for id=1:size(tally,1)
    if strcmpi(tally{id,1},type)==1
        thisid=id;
    end
end

if thisid==0
    tally{idx,1}=type;
else
    type_idx = find(ismember(tally{:,1}, type));
    tally{idx,1}=type;
    tally{idx,2}
end


%%
clear tp;
tp = TypeHierarchy();
b={'F' 'psi' 'KW' 'lbs' 'alarm' 'boolean' 'position'};
F={'MAX_OAT' 'OAT' 'ORH' 'SLCT_PID' 'ART' 'ARS' 'ASO' 'AGN' 'VR' 'MAT' 'SAS' 'SAT' 'RM_SAS' 'RAT' 'SWS' 'SWT' 'CDRWT' 'BLD_1_OAT' 'BLD___ORH'};
lbs = {'SFM' 'HPSTM'};
psi={'HPS' 'VAV__AVG' 'VAV__MIN' 'VAV__MAX' 'VAV' 'RVAV' 'CVP' 'CLV' 'DMP' 'HVP'};
KW={'KWD' 'KWH' 'SKWH' 'KW' 'BLD_1SKWH' 'BLD_1_KWH' 'BLD_1_KWD' 'BLD_2SKWH' 'BLD_2_KWH' 'BLD_2_KWD'};
alarm={'PRALM' 'LOW_RAT' 'LOW_RAT1' 'LOW_RAT2' 'LOW_RAT3' 'LOW_RAT4' 'SMK_ALM' 'SMK_ALM1' 'SMK_ALM2' 'SMK_ALM3' 'SMK_ALM4' 'FLT'};
boolean={'BLD_EVENT' 'BLD___HPS' 'BLD_CURTL' 'BLD_STMON' 'BLD_PRALM'};
position={'VLV' 'A_M' 'L_L'};

tp.addAll('root',b);
tp.addAll('F',F);
tp.addAll('lbs', lbs);
tp.addAll('psi',psi);
tp.addAll('KW',KW);
tp.addAll('alarm', alarm);
tp.addAll('boolean', boolean);
tp.addAll('position',position);
















