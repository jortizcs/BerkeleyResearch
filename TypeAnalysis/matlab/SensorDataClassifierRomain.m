classdef SensorDataClassifierRomain < handle
    %SensorDataClassifier This class uses sensor data to build a
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
        CORY=5
        
        renameAll = '/home/romain/Unison/Projets/jorge/BerkeleyResearch/TypeAnalysis/rnall2';
        
        IMF_bands = ['Hi','Med','Low','Res','Raw'];
        
        tp=0;
        
    end
    
    methods
        function obj=SensorDataClassifierRomain(type, root)

            addpath('/home/romain/Unison/Projets/jorge/BerkeleyResearch/TypeAnalysis/matlab');
            addpath('/home/romain/Unison/Projets/jorge/BerkeleyResearch/TypeAnalysis/matlab/acf');
            
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
            elseif type==obj.CORY
                obj.sourceType=obj.CORY;
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
        
      
        function [accuracy, clusters, DATA] = assess(obj)
            accuracy=0;
            clusters=cell(0,0);
            switch type
                case obj.SODA
                    [accuracy, clusters] = assessSoda(obj);
                case obj.KETI
                    [accuracy, clusters] = assessKETI(obj);
                case obj.SDH
                    [accuracy, clusters] = assessSDH(obj);
                case obj.TODAI
                    [accuracy, clusters] = assessTodai(obj);
            end
        end
        
        % Assess the Soda Hall traces
        % atype is one of the following ['Hist', 'Dpca', 'Raw']
        function [accuracy, pca_acc, clusters] = assessSoda(obj,dict,atype,iband)
            accuracy=0;
            clusters=[];
            if strcmpi(atype,'hist') && obj.sourceType==obj.SODA
                [DATA, GT] = HistDataSoda(obj,dict);
                
            elseif strcmpi(atype,'emd') && obj.sourceType==obj.SODA
                %         IMF_bands = ['Hi','Med','Low','Res','Raw'];

                imfband=5;
                for i=1:length(obj.IMF_bands)
                    if strcmpi(obj.IMF_bands(i),iband)==1
                        imfband=i;
                        break;
                    end
                end
                
                [DATA, GT]=DynDataSoda(obj, dict, 1, imfband);
            else
                fprintf('Valid types, one of: [''Hist'', ''Dpca'', ''Raw'']');
            end
            
            if ~isempty(DATA) && ~isempty(GT)
                mid = floor(size(DATA,1)/2)
                len = length(DATA)
                size(DATA);
                %mdl = ClassificationKNN.fit(DATA(1:mid,:),GT(1:mid));
                idx = knnsearch(DATA(1:mid,:),DATA(:,:),'K',10);
                cnt=0;
                for i=mid+1:len
%                    if strcmp(predict(mdl,DATA(i,:)),GT(i))==1
                if strcmp(GT(idx(i)),GT(i))==1
                       cnt=cnt+1;
                   end
                end

                accuracy=cnt/(len-mid);
                clusters = GT;
                fprintf('\n\nBefore PCA:  accuracy=%f\n',accuracy);
                
                fprintf('Trying w/PCA...');
                % try with PCA
%                 [wcoeff,score,latent,tsquared,explained] = pca(DATA, 'VariableWeights','variance');

                [pc, score,latent] = princomp(DATA);   
                top_k = 1;
                explained = cumsum(latent)./sum(latent);
                while explained(top_k)<0.95
                    top_k=top_k+1;
                end
                %mdl2 = ClassificationKNN.fit(score(1:mid, 1:top_k), GT(1:mid));
                idx2 = knnsearch(score(1:mid,1:top_k),score(:,1:top_k));
                cnt=0;
                for i=mid+1:len
                   %if strcmp(predict(mdl2,score(i,1:top_k)),GT(i))==1
                   if strcmp(GT(idx2(i)),GT(i))==1
                       cnt=cnt+1;
                   end
                end

                pca_acc=cnt/(len-mid);
                fprintf('After PCA:  accuracy=%f\n',pca_acc);

            end
            fprintf('\n\n');
        end
        
        % Assess the KETI traces
        % atype is one of the following ['Hist', 'Dpca', 'Raw']
        function [accuracy, clusters] = assessKETI(obj, atype)
            accuracy=0;
            clusters=[];
            if strcmpi(atype,'hist') && obj.sourceType==obj.KETI
                [DATA, GT] = HistDataKETI(obj);
                mid = floor(size(DATA,1)/5);
                len = size(DATA,1);
                %mdl = ClassificationKNN.fit(DATA(1:mid,:),GT(1:mid));
               %                 idx = knnsearch(DATA(1:mid,:),DATA(:,:),'k',1);
                cnt=0;
                unk=0;
                for i=mid+1:len
                   %search for the hist in the GT that is closer
                   nn = -1;
                   nndiv = 100000;
                   knn = zeros(mid);
                   for j = 1:mid
                    div = KLDiv(DATA(i,:),DATA(j,:));
                    knn(j)=div;
                    if div < nndiv
                        nn = j;
                        nndiv=div;
                    end
                    [tmpX,tmpI]=sort(knn);
                    GT(tmpI(1:3))
                   end
                   if nn == -1
                       unk = unk+1;
                   elseif GT(nn)==GT(i)
                       cnt=cnt+1;
                   end
                end

                accuracy=cnt/(len-mid);
                unkRate=unk/(len-mid)
                clusters = GT;
            else
                fprintf('Valid types, one of: [''Hist'', ''Dpca'', ''Raw'']');
            end
        end
        
        
        % Runs through each KETI trace and reduces it down to
        % a 2-element vector where each element has the height of the
        % first and second highest bars in the histrogram
        % Also returns the ground truth label as follows
        %
        % PIR=3;
        % CO2=1;
        % TEMP=4;
        % HUM=2;
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

                    N=hist(vec,bins);
                    DATA=[DATA; N]; %X(bins) X(bins-1) X(bins-2)];
                    fprintf('\tsize(DATA)=[%d %d]\n\n', size(DATA,1), size(DATA,2));
                end
            end

            fprintf('Done\n\n\n');

        end
        
        % Runs through each UCB Soda Hall trace and reduces it down to
        % a 2-element vector where each element has the height of the
        % first and second highest bars in the histrogram
        % Also returns the ground truth label as follows
        %
        % dict := a cell that has all the labels in it
        function [DATA, GT]=HistDataSoda(obj, dict)
            sz = 5000;
            %base = '/Users/jortiz/Dropbox/dissertation/BerkeleyResearch/TypeAnalysis/data/';
            base = obj.rootdir;
            cd (base);
            folders = dir(base);
            w = 100;    %window for dpca
            DATA=[];

            %ground truth
            GT={};
            idx =1;

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
                        len = 5000;
                        if ~isempty(strfind(files(j,1).name,'M.DAT')) && strcmp(type,'NONE')==0

                            fname = strcat(strcat(thisdir,'/'),files(j,1).name);
                            [stat, res] = system(strcat(['wc -l' ' ' fname]));
                            tok = strtok(res);
                            tlines = str2num(tok);
                            offset = tlines-offset;
                            artcmd = strcat(['ack ''(\d+),\t.*,\t(\d+\.\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);

                            [stat, res]=system(artcmd);
                            artcmd2 = strcat(['ack ''(\\d+),\\t.*,\\t(\\d+\\.\\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);
                            %fprintf(strcat (['artcmd=' artcmd2 '\n']));
                            data = importdata('testres.csv');
                            system('rm -f testres.csv rnall2');

                            % bin and populate
                            if ~isempty(data)
                                GT{idx,1} = type;
                                
                                [N,X]=hist(data(:,2),bins);
                                [~,I]=sort(N);
                                
                                plot(data(:,2));
                                pause
                                %X= sort(X)

                                ac=acf(data(:,2), floor(length(data(:,2))/2));
                                min_w=find(abs(ac(:,1))==min(abs(ac(:,1))));
                                
                                if ~isempty(min_w)
                                    fprintf('\tfile=%s, type=%s, vector=[%d %d %d %d], %d pts, minAC_win=%d\n', ...
                                            files(j,1).name, type, X(bins), X(bins-1), min_w/length(data(:,2)), ac(min_w), size(data,1), min_w);

                                    DATA=[DATA; mean(data(:,2)) std(data(:,2)) min(data(:,2)) max(data(:,2))];
                                end
                                
                                idx = idx+1;
                            end

                        end
                    end
                end


            end

        end
        
        function [DATA, GT]=DynDataSoda(obj, dict, run_emd, imfband)
            
            base = obj.rootdir;
            cd (base);
            sz=5000;
            folders = dir(base);
            w = 100;    %window for dpca
            DATA=[];
            total_pts = 0;

            % number of learning examples per type
            lex=3;
            dpts_per_ex = sz;
            tally = {};
            
            %ground truth
            GT={};
            idx =1;

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
                            artcmd = strcat(['ack ''(\d+),\t.*,\t(\d+\.\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);

                            [stat, res]=system(artcmd);
                            artcmd2 = strcat(['ack ''(\\d+),\\t.*,\\t(\\d+\\.\\d+)$''' ' ' fname ' ' '--output=''$1,$2'' | tail -n ' num2str(offset) ' ' '|head -n' ' ' num2str(len) ' >testres.csv']);
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
                                        newv=min(data(:,1)):15:max(data(:,1));
                                        u=unique(data(:,1));
                                        n=histc(data(:,1),u);
                                        remove_set = find(n>1);
                                        while ~isempty(remove_set)
                                            data(remove_set,:)=[];
                                            fprintf('\t\tDataLength=%d', size(data,1));
                                            newv=min(data(:,1)):15:max(data(:,1));
                                            u=unique(data(:,1));
                                            n=histc(data(:,1),u);
                                            remove_set = find(n>1);
                                        end
                                        % resample the underlying data
                                        data = interp1(data(:,1),data(:,2),newv', 'spline');

                                        % emd shit
                                        [IMFS, REAGG_IMFS] = StripAgg(data,1/15);
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

        
        
         % Assess the SDH traces
        function [accuracy, clusters] = assessSDH(obj)
        end
        
        
        
        % Assess the Todai traces
        function [accuracy, clusters] = assessTodai(obj,atype)
            accuracy=0;
            clusters=[];
            if strcmpi(atype,'hist') && obj.sourceType==obj.TODAI
                [DATA, GT] = HistDataTodai(obj);
                mid = floor(size(DATA,1)/2);
                len = size(DATA,1);
                %mdl = ClassificationKNN.fit(DATA(1:mid,:),GT(1:mid));
               
%                 idx = knnsearch(DATA(1:mid,:),DATA(:,:),'k',1);
                cnt=0;
                unk=0;
                for i=mid+1:len
                   %search for the hist in the GT that is closer
                   nn = -1;
                   nndiv = 100000;
                   for j = 1:mid
                    div = KLDiv(DATA(i,:),DATA(j,:));
                    if div < nndiv
                        nn = j;
                        nndiv=div;
                    end
                   end
                   if nn == -1
                       unk = unk+1;
                   elseif GT(nn)==GT(i)
                       cnt=cnt+1;
                   end
                end

                accuracy=cnt/(len-mid);
                unkRate=unk/(len-mid)
                clusters = GT;
                
            else
                fprintf('Valid types, one of: [''Hist'', ''Dpca'', ''Raw'']');
            end
        end
        
        
        
        % Runs through each Todai trace and reduces it down to
        % a 2-element vector where each element has the height of the
        % first and second highest bars in the histrogram
        % Also returns the ground truth label as follows
        %
% % %         % PIR=3;
% % %         % CO2=1;
% % %         % TEMP=4;
% % %         % HUM=2;
        function [DATA, GT]=HistDataTodai(obj)

            sz = 20000;
            base = obj.rootdir;
            cd (base);
            files = dir([base '/*.dat']);
            DATA=[];
            bins = 10;

            % for ground truth

            EHP=2;
            LIGHT=3;
            POWER=4;
            
%             PIR=3;
%             CO2=1;
%             TEMP=4;
%             HUM=2;
            GT=[];

%             for i=4:size(folder,1)-1
%                 fname= folders(i).name;
%                 files = dir(strcat(base,fname));
%                 thisdir = strcat(strcat(base, fname),'/');
%                 fprintf('Processing %s.\n',thisdir);

                for j=3:size(files,1)
                    fname = files(j).name;

                    % record the type for ground truth
                    type='';
                    if ~isempty(strfind(fname, 'EHP'))
                        GT=[GT; EHP];
                        type='EHP';
                    elseif ~isempty(strfind(fname, 'LIGHT'))
                        GT=[GT; LIGHT];
                        type='LIGHT';
                    elseif ~isempty(strfind(fname, 'POWER'))
                        GT=[GT; POWER];
                        type='POWER';
                    end
                    fprintf('\tGround Truth  [%s]\n',type);

                    fprintf('\tOpening: %s, %d entries\n', strcat(base, fname), sz);
                    fid=fopen(strcat(base, fname));
                    vec = fscanf(fid,'%*s %f',[sz,1]);
                    fclose(fid);

                    N=hist(vec,bins);
                    DATA=[DATA; N]; %X(bins) X(bins-1) X(bins-2)];
                    fprintf('\tsize(DATA)=[%d %d]\n\n', size(DATA,1), size(DATA,2));
                end
%             end

            fprintf('Done\n\n\n');

        end
        
        
        
                % Assess the Cory Hall traces
        function [accuracy, clusters] = assessCory(obj,atype)
            accuracy=0;
            clusters=[];
            if strcmpi(atype,'hist') && obj.sourceType==obj.CORY
                [DATA, GT] = HistDataCory(obj);
                mid = floor(size(DATA,1)/2);
                len = length(DATA);
                %mdl = ClassificationKNN.fit(DATA(1:mid,:),GT(1:mid));
                idx = knnsearch(DATA(1:mid,:),DATA(:,:));
                cnt=0;
                for i=mid+1:len
                   if GT(idx(i))==GT(i)
                       cnt=cnt+1;
                   end
                end

                accuracy=cnt/(len-mid);
                clusters = GT;
            else
                fprintf('Valid types, one of: [''Hist'', ''Dpca'', ''Raw'']');
            end
        end
        
        
        
        % Runs through each Cory Hall trace and reduces it down to
        % a 2-element vector where each element has the height of the
        % first and second highest bars in the histrogram
        % Also returns the ground truth label as follows
        %
% % %         % PIR=3;
% % %         % CO2=1;
% % %         % TEMP=4;
% % %         % HUM=2;
        function [DATA, GT]=HistDataCory(obj)

            sz = 20000;
            base = obj.rootdir;
            cd (base);
            files = dir([base '/*.dat']);
            DATA=[];
            bins = 10;

            % for ground truth
            HVAC=1;
            REC=2;
            ELEV=3;
            LIGHT=4;
            COOL=5;
            
%             PIR=3;
%             CO2=1;
%             TEMP=4;
%             HUM=2;
            GT=[];

%             for i=4:size(folder,1)-1
%                 fname= folders(i).name;
%                 files = dir(strcat(base,fname));
%                 thisdir = strcat(strcat(base, fname),'/');
%                 fprintf('Processing %s.\n',thisdir);

                for j=3:size(files,1)
                    fname = files(j).name;

                    % record the type for ground truth
                    type='';
                    if ~isempty(strfind(fname, 'Cooler'))
                        GT=[GT; COOL];
                        type='Cooler';
                    elseif ~isempty(strfind(fname, 'Hvac'))
                        GT=[GT; HVAC];
                        type='Hvac';
                    elseif ~isempty(strfind(fname, 'Receptacle'))
                        GT=[GT; REC];
                        type='Receptacle';
                    elseif ~isempty(strfind(fname, 'Elevator'))
                        GT=[GT; ELEV];
                        type='Elevator';
                    elseif ~isempty(strfind(fname, 'Light'))
                        GT=[GT; LIGHT];
                        type='Light';
                    end
                    fprintf('\tGround Truth  [%s]\n',type);

                    fprintf('\tOpening: %s, %d entries\n', strcat(base, fname), sz);
                    fid=fopen(strcat(base, fname));
                    vec = fscanf(fid,'%*s %f',[sz,1]);
                    fclose(fid);

                    [N,X]=hist(vec,bins);
                    X= sort(X);
                    DATA=[DATA; mean(vec) std(vec)];
                    fprintf('\tsize(DATA)=[%d %d]\n\n', size(DATA,1), size(DATA,2));
                end
%             end

            fprintf('Done\n\n\n');

        end
    end
    
    
    
    
    
end









