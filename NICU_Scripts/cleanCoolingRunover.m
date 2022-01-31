function cleanCoolingRunover(baseDir)
    hourSegmentsToUse = 2;
    files = dir(fullfile(baseDir, '*.EDF'));
    EEG(1:length(files)) = {eeg_emptyset};
    flatlines=cell(10,1);
    timeBetweenLast = cell(length(files),1);
    if ~exist(fullfile(baseDir,'Output_Cooling_Files'),'dir')
        mkdir(fullfile(baseDir,'Output_Cooling_Files'));
    end
    tic;
    outputDir = 'Output_Cooling_Files';
    delete(fullfile(baseDir,outputDir,'/*.set'));
    delete(fullfile(baseDir,outputDir,'/*.fdt'));
    delete(fullfile(baseDir,'Cut_Images','/*.jpg'));
    fileInfo=table('Size',[length(files), 3],'VariableTypes',{'string','duration','datetime'},'VariableNames',{'Date','Time','Date_Time'});
    fileInfo.Date_Time.Format = 'MM/dd/yyyy H:mm:ss';
    for i=1:length(files)
        fileInfo.Date(i) = datestr(regexp(files(i).name,'\d{2}\-\d{2}\-\d{4}','match'),23);
        fileInfo.Time(i) = duration(regexprep(regexp(files(i).name,'\d{2}\-\d{2}\-\d{2}(?=\.)','match'),{'_','-'},{'',':'}));
        fileInfo.Date_Time(i) = datetime([char(fileInfo.Date(i)) ' ' char(fileInfo.Time(i))],'Format','MM/dd/yyyy H:mm:ss');
 
        EEG{i}=pop_biosig(char(fullfile(baseDir,files(i).name)),'importmex','on');
        EEG{i}=eeg_checkset(EEG{i});

        fileInfo=sortrows(fileInfo,3);
        if i>1
            timeBetweenLast{i} = fileInfo.Date_Time(i)-fileInfo.Date_Time(i-1);
        else
            timeBetweenLast{i} = duration(0,0,0);
        end
        sampRates(i) = EEG{i}.srate;
    end

    %Initialize three 24 hour datasets to be fitted
    if sum(sampRates==EEG{1}.srate) < length(EEG)
        error('One or more of the files do not have identical sampling rate as other files');
    end

    [baseEEG80] = initializeWindows(EEG{1}.srate,EEG{1}.nbchan,EEG{1}.chanlocs);
    
    
    endTime = fileInfo.Date_Time(1) + hours(80);

    for i=1:length(EEG)
        allType={EEG{i}.event(:).type};
        allEdfType={EEG{i}.event(:).edftype};
        allUrEvent={EEG{i}.event(:).urevent};
        allLatency={EEG{i}.event(:).latency};
        if i==1
            pointEnd = EEG{i}.pnts;
            baseEEG80.data(:,1:pointEnd) = EEG{i}.data;
            [baseEEG80.event(1,1:length(EEG{i}.event)).type]=allType{:};
            [baseEEG80.event(1,1:length(EEG{i}.event)).edftype]=allEdfType{:};
            [baseEEG80.event(1,1:length(EEG{i}.event)).urevent]=allUrEvent{:};
            [baseEEG80.event(1,1:length(EEG{i}.event)).latency]=allLatency{:};
            baseEEG80 = eeg_checkset(baseEEG80,'eventconsistency');
        else
            
            pointStart = sum([timeBetweenLast{1:i}]);
            pointStart=hours(pointStart)*60*60*baseEEG80.srate;
            pointEnd = (EEG{i}.pnts-1)+pointStart;
            allLatency = num2cell(cellfun(@(x) x+pointStart,allLatency));
            baseEEG80.data(:,pointStart:pointEnd) = EEG{i}.data;
            if length(baseEEG80.event(:))>0
                startIndex = length(baseEEG80.event(:))+1;
            else
                startIndex=1;
            end
            [baseEEG80.event(1,startIndex:startIndex+length(EEG{i}.event)-1).type]=allType{:};
            [baseEEG80.event(1,startIndex:startIndex+length(EEG{i}.event)-1).edftype]=allEdfType{:};
            [baseEEG80.event(1,startIndex:startIndex+length(EEG{i}.event)-1).urevent]=allUrEvent{:};
            [baseEEG80.event(1,startIndex:startIndex+length(EEG{i}.event)-1).latency]=allLatency{:};
            baseEEG80 = eeg_checkset(baseEEG80,'eventconsistency');
        end

    end
    toc
    EEG=[];
    %Binning and flatline
    for i=1:10
        switch i
            case 1
                baseEEG8=eeg_emptyset;
                baseEEG8 = pop_select(baseEEG80,'point',[1 8*60*60*baseEEG80.srate],'channel',1:baseEEG80.nbchan);
                [flatlines{i},baseEEG8] = findFlatlines(baseEEG8,baseDir,1,hourSegmentsToUse);
                baseEEG8.etc.processStep = 'Cooling';
                baseEEG8.etc.hoursBin=1;
                baseEEG8.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG8.srate)])>0))
                    pop_saveset(baseEEG8,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG8.etc.processStep])));
                else
                    baseEEG8=[];
                end
            case 2
                baseEEG16=eeg_emptyset;
                baseEEG16 = pop_select(baseEEG80,'point',[8*60*60*baseEEG80.srate+1 8*2*60*60*baseEEG80.srate],'channel',1:baseEEG80.nbchan);
                [flatlines{i},baseEEG16] = findFlatlines(baseEEG16,baseDir,2,hourSegmentsToUse);
                baseEEG16.etc.processStep = 'Cooling';
                baseEEG16.etc.hoursBin=2;
                baseEEG16.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG16.srate)])>0))
                    pop_saveset(baseEEG16,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG16.etc.processStep])));
                else
                    baseEEG16=[];
                end
            case 3
                baseEEG24=eeg_emptyset;
                baseEEG24 = pop_select(baseEEG80,'point',[8*2*60*60*baseEEG80.srate+1 8*3*60*60*baseEEG80.srate],'channel',1:baseEEG80.nbchan);
                [flatlines{i},baseEEG24] = findFlatlines(baseEEG24,baseDir,3,hourSegmentsToUse);
                baseEEG24.etc.processStep = 'Cooling';
                baseEEG24.etc.hoursBin=3;
                baseEEG24.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG24.srate)])>0))
                    pop_saveset(baseEEG24,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG24.etc.processStep])));
                else
                    baseEEG24=[];
                end
            case 4
                baseEEG32=eeg_emptyset;
                baseEEG32 = pop_select(baseEEG80,'point',[8*3*60*60*baseEEG80.srate+1 8*4*60*60*baseEEG80.srate],'channel',1:baseEEG80.nbchan);
                [flatlines{i},baseEEG32] = findFlatlines(baseEEG32,baseDir,4,hourSegmentsToUse);
                baseEEG32.etc.processStep = 'Cooling';
                baseEEG32.etc.hoursBin=4;
                baseEEG32.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG32.srate)])>0))
                    pop_saveset(baseEEG32,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG32.etc.processStep])));
                else
                    baseEEG32=[];
                end
            case 5
                baseEEG40=eeg_emptyset;
                baseEEG40 = pop_select(baseEEG80,'point',[8*4*60*60*baseEEG80.srate+1 8*5*60*60*baseEEG80.srate],'channel',1:baseEEG80.nbchan);
                [flatlines{i},baseEEG40] = findFlatlines(baseEEG40,baseDir,5,hourSegmentsToUse);
                baseEEG40.etc.processStep = 'Cooling';
                baseEEG40.etc.hoursBin=5;
                baseEEG40.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG40.srate)])>0))
                    pop_saveset(baseEEG40,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG40.etc.processStep])));
                else
                    baseEEG40=[];
                end
            case 6
                baseEEG48 = eeg_emptyset;
                baseEEG48 = pop_select(baseEEG80,'point',[8*5*60*60*baseEEG80.srate+1 8*6*60*60*baseEEG80.srate],'channel',1:baseEEG80.nbchan);
                [flatlines{i},baseEEG48] = findFlatlines(baseEEG48,baseDir,6,hourSegmentsToUse);
                baseEEG48.etc.processStep = 'Cooling';
                baseEEG48.etc.hoursBin=6;
                baseEEG48.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG48.srate)])>0))
                    pop_saveset(baseEEG48,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG48.etc.processStep])));
                else
                    baseEEG48=[];
                end
            case 7
                baseEEG56 = eeg_emptyset;
                baseEEG56 = pop_select(baseEEG80,'point',[8*6*60*60*baseEEG80.srate+1 8*7*60*60*baseEEG80.srate],'channel',1:baseEEG80.nbchan);
                [flatlines{i},baseEEG56] = findFlatlines(baseEEG56,baseDir,7,hourSegmentsToUse);
                baseEEG56.etc.processStep = 'Cooling';
                baseEEG56.etc.hoursBin=7;
                baseEEG56.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG56.srate)])>0))
                    pop_saveset(baseEEG56,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG56.etc.processStep])));
                else
                    baseEEG56=[];
                end
            case 8
                baseEEG64 = eeg_emptyset;
                baseEEG64 = pop_select(baseEEG80,'point',[8*7*60*60*baseEEG80.srate+1 8*8*60*60*baseEEG80.srate],'channel',1:baseEEG80.nbchan);
                [flatlines{i},baseEEG64] = findFlatlines(baseEEG64,baseDir,8,hourSegmentsToUse);
                baseEEG64.etc.processStep = 'Cooling';
                baseEEG64.etc.hoursBin=8;
                baseEEG64.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG64.srate)])>0))
                    pop_saveset(baseEEG64,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG64.etc.processStep])));
                else
                    baseEEG64=[];
                end
            case 9
                baseEEG72 = pop_select(baseEEG80,'point',[8*8*60*60*baseEEG80.srate+1 8*9*60*60*baseEEG80.srate],'channel',1:baseEEG80.nbchan);
                [flatlines{i},baseEEG72] = findFlatlines(baseEEG72,baseDir,9,hourSegmentsToUse);
                baseEEG72.etc.processStep = 'Cooling';
                baseEEG72.etc.hoursBin=9;
                baseEEG72.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG72.srate)])>0))
                    pop_saveset(baseEEG72,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG72.etc.processStep])));
                else
                    baseEEG72=[];
                end
                
            case 10
                 baseEEG80 = pop_select(baseEEG80,'point',[9*8*60*60*baseEEG80.srate+1 8*10*60*60*baseEEG80.srate],'channel',1:baseEEG80.nbchan);
                [flatlines{i},baseEEG80] = findFlatlines(baseEEG80,baseDir,10,hourSegmentsToUse);
                baseEEG80.etc.processStep = 'Cooling';
                baseEEG80.etc.hoursBin=10;
                baseEEG80.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG80.srate)])>0))
                    pop_saveset(baseEEG80,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG80.etc.processStep])));
                else
                    baseEEG80=[];
                end
        end
    end

    function [baseEEG80]=initializeWindows(srate,nbchan,chanlocs)
        baseEEG80 = eeg_emptyset();

        baseEEG80.srate = srate;
        baseEEG80.nbchan = nbchan;
        baseEEG80.chanlocs = chanlocs;
        baseEEG80.data(1:baseEEG80.nbchan,1:(80*60*60)*baseEEG80.srate)=0;
        baseEEG80.event = struct('type', {}, 'edftype',{},'urevent',{},'latency',{});
        baseEEG80 = eeg_checkset(baseEEG80);
        

    end
end