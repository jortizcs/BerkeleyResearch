classdef SensorDataAnalyzer < handle
    %SensorDataAnalyzer This class uses sensor data to build a
    %classifier.
    %   Classification is performed using several approaches/options.  We
    %   inlcude DPCA, DPCA+EMD, and several others.
    
    properties (SetAccess = public)
        sourceType=0
        rootdir=zeros(0,0)
        
        NONE=0
        SODA=1
        KETI=2
        SDH=3
        TODAI=4
        
        renameAll = '/Users/jortiz/data/rnall2';
        
        IMF_bands = ['Hi','Med','Low','Res','Raw'];
        
        tp=0;
        filestats = cell(0,4);
    end
    
    methods
        function obj=SensorDataAnalyzer(type, root)

            addpath('/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/matlab');
            addpath('/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/matlab/acf');
            
            % for hierarchical clustering, this sets up the label hiearchy
            obj.tp = TypeHierarchy();
            b={'F' 'psi' 'KW' 'lbs' 'alarm' 'boolean' 'position'};
            F={'MAX_OAT' 'OAT' 'ORH' 'SLCT_PID' 'ART' 'ARS' 'ASO' 'AGN' 'VR' 'MAT' 'SAS' 'SAT' 'RM_SAS' 'RAT' 'SWS' 'SWT' 'CDRWT' 'BLD_1_OAT' 'BLD___ORH'};
            lbs = {'SFM' 'HPSTM'};
            psi={'HPS' 'VAV__AVG' 'VAV__MIN' 'VAV__MAX' 'VAV' 'RVAV' 'CVP' 'CLV' 'DMP' 'HVP'};
            KW={'KWD' 'KWH' 'SKWH' 'KW' 'BLD_1SKWH' 'BLD_1_KWH' 'BLD_1_KWD' 'BLD_2SKWH' 'BLD_2_KWH' 'BLD_2_KWD'};
            alarm={'PRALM' 'LOW_RAT' 'LOW_RAT1' 'LOW_RAT2' 'LOW_RAT3' 'LOW_RAT4' 'SMK_ALM' 'SMK_ALM1' 'SMK_ALM2' 'SMK_ALM3' 'SMK_ALM4' 'FLT'};
            boolean={'BLD_EVENT' 'BLD___HPS' 'BLD_CURTL' 'BLD_STMON' 'BLD_PRALM'};
            position={'VLV' 'A_M' 'L_L'};

            obj.tp.addAll('root',b);
            obj.tp.addAll('F',F);
            obj.tp.addAll('lbs', lbs);
            obj.tp.addAll('psi',psi);
            obj.tp.addAll('KW',KW);
            obj.tp.addAll('alarm', alarm);
            obj.tp.addAll('boolean', boolean);
            obj.tp.addAll('position',position);

            %type setting
            if type==obj.NONE
                obj.sourceType = obj.NONE;
            elseif type==obj.SODA
                obj.sourceType=obj.SODA;
            elseif type==obj.KETI
                obj.sourceType=obj.KETI;
            elseif type==obj.TODAI
                obj.sourceType=obj.TODAI;
            end
           
            if root(1,length(root))~='/'
                root = strcat(root,'/');
            end
            obj.rootdir=root;
        end
        
        function sourceType=setSourceType(obj,type)
            if type==obj.NONE
                obj.sourceType = obj.NONE;
            elseif type==obj.SODA
                obj.sourceType=obj.SODA;
            elseif type==obj.KETI
                obj.sourceType=obj.KETI;
            elseif type==obj.TODAI
                obj.sourceType=obj.TODAI;
            end        
            sourceType=obj.sourceType;
        end
        
        function r=setRootDirectory(obj,r)
            obj.rootdir=r;
        end
        
        % Assess the KETI traces
        % atype is one of the following ['Hist', 'Dpca', 'Raw']
        function [accuracy, clusters] = assessKETI(obj, atype)
            accuracy=0;
            clusters=[];
            if strcmpi(atype,'hist') && obj.sourceType==obj.KETI
                [DATA, GT] = HistDataKETI(obj);
                mid = floor(size(DATA,1)/2);
                len = length(DATA);
                mdl = ClassificationKNN.fit(DATA(1:mid,:),GT(1:mid));

                cnt=0;
                for i=mid+1:len
                   if predict(mdl,DATA(i,:))==GT(i)
                       cnt=cnt+1;
                   end
                end

                accuracy=cnt/(len-mid);
                clusters = GT;
            else
                fprintf('Valid types, one of: [''Hist'', ''Dpca'', ''Raw'']');
            end
        end
        
        function [accuracy] = assessSoda(obj, dict)
            %[DATA, MVARDATA, GT, GTMVAR, stats]=obj.DynDataSoda(dict, 1, 5);
            [DATA, MVARDATA, GT, GTMVAR, stats]=obj.Test_DynDataSoda(dict, 1, 5);
            
            fprintf('Running KNN\n');
            
            if ~isempty(DATA) && ~isempty(GT)
                mid = floor(size(DATA,1)/2)
                len = length(DATA)
                size(DATA)
                size(GT)
                mdl = ClassificationKNN.fit(DATA(1:mid,:),GT(1:mid));
                
                cnt=0;
                for i=mid+1:len
                   if strcmp(predict(mdl,DATA(i,:)),GT(i))==1
                       cnt=cnt+1;
                   end
                end

                accuracy=cnt/(len-mid);
                clusters = GT;
                fprintf('\n\nBefore PCA:  accuracy=%f\n',accuracy);
                
                fprintf('Trying w/PCA...');
                % try with PCA
                [wcoeff,score,latent,tsquared,explained] = pca(DATA, 'VariableWeights','variance');
                top_k = 1;
                while sum(explained(1,1:top_k))<0.95
                    top_k=top_k+1;
                end
                mdl2 = ClassificationKNN.fit(score(1:mid, 1:top_k), GT(1:mid));
                cnt=0;
                for i=mid+1:len
                   if strcmp(predict(mdl2,score(i,1:top_k)),GT(i))==1
                       cnt=cnt+1;
                   end
                end

                pca_acc=cnt/(len-mid);
                fprintf('After PCA:  accuracy=%f\n',pca_acc);

            end
            
            
            fprintf('\n\nAnalyzing mean & standard deviation approach...\n');
            if ~isempty(MVARDATA) && ~isempty(GTMVAR)
                mid = floor(size(MVARDATA,1)/2);
                len = length(MVARDATA);
                mdl = ClassificationKNN.fit(MVARDATA(1:mid,:),GTMVAR(1:mid));
                
                cnt=0;
                for i=mid+1:len
                   if strcmp(predict(mdl,MVARDATA(i,:)),GTMVAR(i))==1
                       cnt=cnt+1;
                   end
                end

                accuracy=cnt/(len-mid);
                clusters = GTMVAR;
                fprintf('\n\nBefore PCA:  accuracy=%f\n',accuracy);
                
                fprintf('Trying w/PCA...');
                % try with PCA
                [wcoeff,score,latent,tsquared,explained] = pca(MVARDATA, 'VariableWeights','variance');
                top_k = 1;
                while sum(explained(1,1:top_k))<0.95
                    top_k=top_k+1;
                end
                mdl2 = ClassificationKNN.fit(score(1:mid, 1:top_k), GTMVAR(1:mid));
                cnt=0;
                for i=mid+1:len
                   if strcmp(predict(mdl2,score(i,1:top_k)),GTMVAR(i))==1
                       cnt=cnt+1;
                   end
                end

                pca_acc=cnt/(len-mid);
                fprintf('After PCA:  accuracy=%f\n',pca_acc);

            end
            fprintf('\n\n');
        end
        
        
        % Runs through each KETI trace and reduces it down to
        % a 2-element vector where each element has the height of the
        % first and second highest bars in the histrogram
        % Also returns the ground truth label as follows
        %
        function [DATA, GT]=HistDataKETI(obj)

            sz = 20000;
            %base = '~/StreamFS/StreamFS.classifier/data/KETI/';
            base = obj.rootdir;
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
0
                    [N,X]=hist(vec,bins);
                    X= sort(X);
                    DATA=[DATA; X(bins) X(bins-1) X(bins-2)];
                    fprintf('\tsize(DATA)=[%d %d]\n\n', size(DATA,1), size(DATA,2));
                end
            end

            fprintf('Done\n\n\n');

        end
        
        
        
        
        
        
        
        
        
        function [DATA, MVARDATA, GT, GTMVAR, stats]=DynDataSoda(obj, dict, run_emd, imfband)
            
            base = obj.rootdir;
            cd (base);
            sz=5000;
            folders = dir(base);
            w = 100;    %window for dpca
            DATA=[];
            total_pts = 0;
            stats = cell(0,0);
            sidx = 1;
            MVARDATA = []; %trace mean and std dev.
            
            % number of learning examples per type
            lex=3;
            dpts_per_ex = sz;
            tally = {};
            
            %ground truth
            GT={};
            idx =1;
            GTMVAR={};

            % bins
            bins=10;

            % figure;
            for i=4:size(folders,1)-1
                fname= folders(i).name;
                thisdir = strcat(strcat(base, fname));
                if isdir(thisdir)        
                    system(strcat(['cp' ' ' obj.renameAll ' ' thisdir]));
                    cd(thisdir);
                    system('source rnall2');
                    fprintf('Processing %s\n',thisdir);

                    % figure out the type
                    type = '';
                    for i=1:length(dict)
                        if ~isempty(strfind(fname, dict{1,i})) && length(dict{1,i})>length(type)
                            type=dict{1,i};
                        end
                    end
                    
                    if strcmp(type,'')==1
                        type='NONE';
                    end

                    % start parsing
                    files = dir(strcat(strcat(base,fname),'/'));
                    for j=3:size(files,1)
                        offset = 0;
                        len = sz;
                        if ~isempty(strfind(files(j,1).name,'M.DAT')) && strcmp(type,'NONE')==0

                            fprintf('\tProcessing: %s\n', files(j,1).name);
                            fname = strcat(strcat(thisdir,'/'),files(j,1).name);
                            [stat, res] = system(strcat(['wc -l' ' ' fname]));
                            tok = strtok(res);
                            tlines = str2num(tok);
                            offset = tlines-offset;
                            artcmd = strcat(['~/bin/ack ''(\d+),\t.*,\t(\d+\.\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);

                            [stat, res]=system(artcmd);
                            artcmd2 = strcat(['~/bin/ack ''(\\d+),\\t.*,\\t(\\d+\\.\\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);
                            %fprintf(strcat (['artcmd=' artcmd2 '\n']));
                            data = importdata('testres.csv');
                            system('rm -f testres.csv rnall2');

                            % bin and populate
                            if ~isempty(data)
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % keep track of number learning examples%%%
                                do_not_include =0;
                                thisid =0;
                                for id=1:size(tally,1)
                                    if strcmpi(tally{id,1},type)==1
                                        thisid=id;
                                    end
                                end

                                if thisid==0
                                    newid=size(tally,1)+1;
                                    tally{newid,1}=type;
                                    tally{newid,2}=length(data);
                                else
                                    if tally{thisid,2}>dpts_per_ex*lex
                                        do_not_include = 1;
                                    else
                                        tally{thisid,2}=tally{thisid,2}+length(data);
                                    end
                                end
                                tally
                                if size(tally,1)==size(dict,2) && total_pts>=size(dict,2)*dpts_per_ex*lex
                                    fprintf('Done learning\n');
                                    return;
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                % include and analyze it only if not seen
                                % before
                                if do_not_include==0
                                    total_pts = total_pts + length(data);
                                    if run_emd==1
                                        fprintf('\tRunning EMD on %s...\n',files(j,1).name);

                                        %removing deplicates
                                        fprintf('\t\tDataLength=%d\n', size(data,1));
                                        newv=min(data(:,1)):60:max(data(:,1));
                                        u=unique(data(:,1));
                                        n=histc(data(:,1),u);
                                        remove_set = find(n>1);
                                        while ~isempty(remove_set)
                                            data(remove_set,:)=[];
                                            fprintf('\t\tDataLength=%d', size(data,1));
                                            newv=min(data(:,1)):60:max(data(:,1));
                                            u=unique(data(:,1));
                                            n=histc(data(:,1),u);
                                            remove_set = find(n>1);
                                        end
                                        
                                        iratio = length(data);
                                        % resample the underlying data
                                        data = interp1(data(:,1),data(:,2),newv', 'spline');
                                        iratio= length(data)/iratio
                                        stats{sidx,1}=files(j,1).name;
                                        stats{sidx,2} = {'interpolation ratio', iratio};
                                        sidx = sidx +1;
                                        
                                        mn = mean(data(:,1))
                                        sd = std(data(:,1))
                                        type
                                        MVARDATA = [MVARDATA; mn, sd];
                                        GTMVAR{length(GTMVAR)+1,1} = type;

                                        % emd shit
                                        [IMFS, REAGG_IMFS] = StripAgg(data,1/60);
                                        size(REAGG_IMFS)
                                        imfband
                                        data = REAGG_IMFS(imfband,:);
                                        fprintf('\t...done\n');
                                    end

                                    chunks = floor(length(data)/w);
                                    fprintf('\t\tData not empty, chunks=%d\n', chunks);

%                                     ac=acf(data(:,2), floor(length(data(:,2))/2));
%                                     min_w=find(abs(ac(:,1))==min(abs(ac(:,1))));

                                    for i=1:chunks
                                        %size(data);
                                        DATA=[DATA; data(1,1:w)];
                                        GT{idx,1}=type;
                                        idx = idx +1;
                                        data = data(1,w+1:length(data));
                                    end
%                                     fprintf('\t\t\tfile=%s, type=%s, vector=[%d %d %d ...], %d pts, minAC_win=%d\n', ...
%                                             files(j,1).name, type, data(1,1), data(1,2), data(1,3), size(data,2), min_w);
                                end
                                
                            end
%                         else
%                             fprintf('\tInvalid: %s\n', files(j,1).name);
                        end
                    end
                end

            end
        end
        
        
         function [DATA, MVARDATA, GT, GTMVAR, stats]=Test_DynDataSoda(obj, dict, run_emd, imfband)
            
            base = obj.rootdir;
            cd (base);
            sz=20000;
            folders = dir(base);
            w = 100;    %window for dpca
            DATA=[];
            total_pts = 0;
            stats = cell(0,0);
            sidx = 1;
            MVARDATA = []; %trace mean and std dev.
            
            % number of learning examples per type
            lex=3;
            dpts_per_ex = sz;
            tally = {};
            
            %ground truth
            GT={};
            idx =1;
            GTMVAR={};

            % bins
            bins=10;

            % figure;
            for i=4:size(folders,1)-1
                fname= folders(i).name;
                thisdir = strcat(strcat(base, fname));
                if isdir(thisdir)        
                    system(strcat(['cp' ' ' obj.renameAll ' ' thisdir]));
                    cd(thisdir);
                    system('source rnall2');
                    fprintf('Processing %s\n',thisdir);

                    % figure out the type
                    type = '';
                    for i=1:length(dict)
                        if ~isempty(strfind(fname, dict{1,i})) && length(dict{1,i})>length(type)
                            type=dict{1,i};
                        end
                    end
                    
                    if strcmp(type,'')==1
                        type='NONE';
                    end

                    % start parsing
                    files = dir(strcat(strcat(base,fname),'/'));
                    for j=3:size(files,1)
                        offset = 0;
                        len = sz;
                        if ~isempty(strfind(files(j,1).name,'M.DAT')) && strcmp(type,'NONE')==0

                            fprintf('\tProcessing: %s\n', files(j,1).name);
                            fname = strcat(strcat(thisdir,'/'),files(j,1).name);
                            [stat, res] = system(strcat(['wc -l' ' ' fname]));
                            tok = strtok(res);
                            tlines = str2num(tok);
                            offset = tlines-offset;
                            artcmd = strcat(['~/bin/ack ''(\d+),\t.*,\t(\d+\.\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);

                            [stat, res]=system(artcmd);
                            artcmd2 = strcat(['~/bin/ack ''(\\d+),\\t.*,\\t(\\d+\\.\\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);
                            %fprintf(strcat (['artcmd=' artcmd2 '\n']));
                            data = importdata('testres.csv');
                            system('rm -f testres.csv rnall2');

                            % bin and populate
                            if ~isempty(data)
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % keep track of number learning examples%%%
                                do_not_include =0;
                                thisid =0;
                                for id=1:size(tally,1)
                                    if strcmpi(tally{id,1},type)==1
                                        thisid=id;
                                    end
                                end

                                if thisid==0
                                    newid=size(tally,1)+1;
                                    tally{newid,1}=type;
                                    tally{newid,2}=length(data);
                                else
                                    if tally{thisid,2}>dpts_per_ex*lex
                                        do_not_include = 1;
                                    else
                                        tally{thisid,2}=tally{thisid,2}+length(data);
                                    end
                                end
                                tally
                                if size(tally,1)==size(dict,2) && total_pts>=size(dict,2)*dpts_per_ex*lex
                                    fprintf('Done learning\n');
                                    return;
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                % include and analyze it only if not seen
                                % before
                                if do_not_include==0
                                    total_pts = total_pts + length(data);
                                    if run_emd==1
                                        fprintf('\tRunning EMD on %s...\n',files(j,1).name);

                                        %removing deplicates
                                        fprintf('\t\tDataLength=%d\n', size(data,1));
                                        newv=min(data(:,1)):60:max(data(:,1));
                                        u=unique(data(:,1));
                                        n=histc(data(:,1),u);
                                        remove_set = find(n>1);
                                        while ~isempty(remove_set)
                                            data(remove_set,:)=[];
                                            fprintf('\t\tDataLength=%d', size(data,1));
                                            newv=min(data(:,1)):60:max(data(:,1));
                                            u=unique(data(:,1));
                                            n=histc(data(:,1),u);
                                            remove_set = find(n>1);
                                        end
                                        
                                        obj.filestats{length(obj.filestats)+1,1}=files(j,1).name;
                                        obj.filestats{length(obj.filestats),2} = mean(data(:,2));
                                        obj.filestats{length(obj.filestats),3} = std(data(:,2));
                                        obj.filestats{length(obj.filestats),4} = mode(data(:,2));
                                        
                                        iratio = length(data);
                                        % resample the underlying data
                                        data = interp1(data(:,1),data(:,2),newv', 'spline');
                                        iratio= length(data)/iratio
                                        stats{sidx,1}=files(j,1).name;
                                        stats{sidx,2} = {'interpolation ratio', iratio};
                                        sidx = sidx +1;
                                                                                
                                        mn = mean(data(:,1))
                                        sd = std(data(:,1))
                                        type
                                        MVARDATA = [MVARDATA; mn, sd];
                                        GTMVAR{length(GTMVAR)+1,1} = type;

                                        % emd shit
                                        %[IMFS, REAGG_IMFS] = StripAgg(data,1/60);
%                                         figure;
%                                         subplot(2,1,1);
%                                         plot(REAGG_IMFS(imfband,:));
%                                         subplot(2,1,2);
                                        f=1/60; 

                                        order = 3; % order of the filter is 3

                                        fnorm1 = [1/(20*60)]/(f/2); % for bandpass, here 1 and 3 are the lower and upper cutoff respectively
                                        % [b,a] = butter(n,Wn) % n = order of the filter, wn is the normalized
                                        % cutoff freq
                                        [b1,a1] = butter(order,fnorm1,'high'); 
                                        data_1 = filtfilt(b1,a1,data); % band pass filter
                                        
                                        fnorm2 = [1/(1440*60) 1/(20*60)]/(f/2); % for bandpass
                                        [b2,a2] = butter(order,fnorm2);
                                        data_2 = filtfilt(b2,a2,data); % band pass filter

                                        fnorm3 = [1/(2160*60)]/(f/2); % for bandpass
                                        [b3,a3] = butter(order,fnorm3,'low');
                                        data_3 = filtfilt(b3,a3,data); % band pass filter
                                        %data_4 = data - (data_1+data_2+data_3); % residual

                                        %plot(data_3);
                                        %return;
                                        fprintf('Residual extracted\n');
                                        
%                                         size(REAGG_IMFS)
%                                         imfband
%                                         data = REAGG_IMFS(imfband,:);
                                        data = data_3;
                                        fprintf('\t...done\n');
                                    end

                                    chunks = floor(length(data)/w);
                                    fprintf('\t\tData not empty, chunks=%d\n', chunks);

%                                     ac=acf(data(:,2), floor(length(data(:,2))/2));
%                                     min_w=find(abs(ac(:,1))==min(abs(ac(:,1))));

                                    for i=1:chunks
                                        %size(data);
                                        %DATA=[DATA; data(1,1:w)];
                                        DATA=[DATA; data(1:w)'];
                                        GT{idx,1}=type;
                                        idx = idx +1;
                                        %data = data(1,w+1:length(data));
                                        data = data(w+1:length(data));
                                    end
                                    size(DATA)

%                                     fprintf('\t\t\tfile=%s, type=%s, vector=[%d %d %d ...], %d pts, minAC_win=%d\n', ...
%                                             files(j,1).name, type, data(1,1), data(1,2), data(1,3), size(data,2), min_w);
                                end
                                
                            end
%                         else
%                             fprintf('\tInvalid: %s\n', files(j,1).name);
                        end
                    end
                end

            end
         end
        
         
        function [base]=Test2_DynDataSoda(obj, dict)
            
            base = obj.rootdir;
            cd (base);
            sz=20000;
            folders = dir(base);
            
            % number of learning examples per type
            lex=3;
            dpts_per_ex = sz;
            tally = {};

            % bins
            bins=10;

            % figure;
            for i=4:size(folders,1)-1
                fname= folders(i).name;
                thisdir = strcat(strcat(base, fname));
                if isdir(thisdir)        
                    system(strcat(['cp' ' ' obj.renameAll ' ' thisdir]));
                    cd(thisdir);
                    system('source rnall2');
                    fprintf('Processing %s\n',thisdir);

                    % figure out the type
                    type = '';
                    for i=1:length(dict)
                        if ~isempty(strfind(fname, dict{1,i})) && length(dict{1,i})>length(type)
                            type=dict{1,i};
                        end
                    end
                    
                    if strcmp(type,'')==1
                        type='NONE';
                    end

                    % start parsing
                    files = dir(strcat(strcat(base,fname),'/'));
                    for j=3:size(files,1)
                        offset = 0;
                        len = sz;
                        if ~isempty(strfind(files(j,1).name,'M.DAT')) && strcmp(type,'NONE')==0

                            fprintf('\tProcessing: %s\n', files(j,1).name);
                            fname = strcat(strcat(thisdir,'/'),files(j,1).name);
                            [stat, res] = system(strcat(['wc -l' ' ' fname]));
                            tok = strtok(res);
                            tlines = str2num(tok);
                            offset = tlines-offset;
                            artcmd = strcat(['~/bin/ack ''(\d+),\t.*,\t(\d+\.\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);

                            [stat, res]=system(artcmd);
                            artcmd2 = strcat(['~/bin/ack ''(\\d+),\\t.*,\\t(\\d+\\.\\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);
                            %fprintf(strcat (['artcmd=' artcmd2 '\n']));
                            data = importdata('testres.csv');
                            system('rm -f testres.csv rnall2');

                            % bin and populate
                            if ~isempty(data)
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % keep track of number learning examples%%%
                                do_not_include =0;
                                thisid =0;
                                for id=1:size(tally,1)
                                    if strcmpi(tally{id,1},type)==1
                                        thisid=id;
                                    end
                                end

                                if thisid==0
                                    newid=size(tally,1)+1;
                                    tally{newid,1}=type;
                                    tally{newid,2}=length(data);
                                else
                                    if tally{thisid,2}>dpts_per_ex*lex
                                        do_not_include = 1;
                                    else
                                        tally{thisid,2}=tally{thisid,2}+length(data);
                                    end
                                end
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                % include and analyze it only if not seen
                                % before

                                %removing deplicates
                                fprintf('\t\tDataLength=%d\n', size(data,1));
                                newv=min(data(:,1)):60:max(data(:,1));
                                u=unique(data(:,1));
                                n=histc(data(:,1),u);
                                remove_set = find(n>1);
                                while ~isempty(remove_set)
                                    data(remove_set,:)=[];
                                    fprintf('\t\tDataLength=%d\n', size(data,1));
                                    newv=min(data(:,1)):60:max(data(:,1));
                                    u=unique(data(:,1));
                                    n=histc(data(:,1),u);
                                    remove_set = find(n>1);
                                end

                                mn = mean(data(:,2));
                                sd = std(data(:,2));
                                md = mode(data(:,2));
                                if ~isempty(mn) && ~isempty(sd) && ~isempty(md)
                                    obj.filestats{length(obj.filestats)+1,1}=files(j,1).name;
                                    obj.filestats{length(obj.filestats),2} = mn;
                                    obj.filestats{length(obj.filestats),3} = sd;
                                    obj.filestats{length(obj.filestats),4} = md;

                                    fprintf(strcat(['[' files(j,1).name ', ' num2str(mn) ', ' num2str(sd) ', ' num2str(md) ']\n']));
                                end
                                
                            end
                        end
                    end

                end
            end
         
        end
    
    
    
    
         
        
        
    
    end
         
         
         
         
         
         
         
         
         
    
end









