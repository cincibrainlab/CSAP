function cleanRewarming(baseDir)
    files = dir(fullfile(baseDir, '*.EDF'));
    startInfo = cell(length(files),3);
    EEG(1:length(files)) = {eeg_emptyset};
    flatlines=cell(6,1);
    timeBetweenLast = cell(length(files),1);
    if ~exist(fullfile(baseDir,'Output_Rewarming_Files'),'dir')
        mkdir(fullfile(baseDir,'Output_Rewarming_Files'));
    end
    outputDir = 'Output_Rewarming_Files';
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
    
    if sum(sampRates==EEG{1}.srate) < length(EEG)
        error('One or more of the files do not have identical sampling rate as other files');
    end
    
    %Initialize three 24 hour datasets to be fitted
    [baseEEG48] = initializeWindows(EEG{1}.srate,EEG{1}.nbchan,EEG{1}.chanlocs);
    endTime = fileInfo.Date_Time(1) + hours(48);
    for i=1:length(EEG)
        allType={EEG{i}.event(:).type};
        allEdfType={EEG{i}.event(:).edftype};
        allUrEvent={EEG{i}.event(:).urevent};
        allLatency={EEG{i}.event(:).latency};
        if i==1
            pointEnd = EEG{i}.pnts;
            baseEEG48.data(:,1:pointEnd) = EEG{i}.data;
            [baseEEG48.event(1,1:length(EEG{i}.event)).type]=allType{:};
            [baseEEG48.event(1,1:length(EEG{i}.event)).edftype]=allEdfType{:};
            [baseEEG48.event(1,1:length(EEG{i}.event)).urevent]=allUrEvent{:};
            [baseEEG48.event(1,1:length(EEG{i}.event)).latency]=allLatency{:};
            baseEEG48 = eeg_checkset(baseEEG48,'eventconsistency');
        else
            
            pointStart = sum([timeBetweenLast{1:i}]);
            pointStart=hours(pointStart)*60*60*baseEEG48.srate;
            pointEnd = (EEG{i}.pnts-1)+pointStart;
            allLatency = num2cell(cellfun(@(x) x+pointStart,allLatency));
            baseEEG48.data(:,pointStart:pointEnd) = EEG{i}.data;
            if length(baseEEG48.event(:))>0
                startIndex = length(baseEEG48.event(:))+1;
            else
                startIndex=1;
            end
            [baseEEG48.event(1,startIndex:startIndex+length(EEG{i}.event)-1).type]=allType{:};
            [baseEEG48.event(1,startIndex:startIndex+length(EEG{i}.event)-1).edftype]=allEdfType{:};
            [baseEEG48.event(1,startIndex:startIndex+length(EEG{i}.event)-1).urevent]=allUrEvent{:};
            [baseEEG48.event(1,startIndex:startIndex+length(EEG{i}.event)-1).latency]=allLatency{:};
            baseEEG48 = eeg_checkset(baseEEG48, 'eventconsistency');
        end
    end
% end
    EEG=[];
    %Binning and flatline
    for i=1:6
        switch i
            case 1
                baseEEG8=eeg_emptyset;
                baseEEG8 = pop_select(baseEEG48,'point',[1 8*60*60*baseEEG48.srate],'channel',1:baseEEG48.nbchan);
                [flatlines{i},baseEEG8] = findFlatlines(baseEEG8,baseDir,1,2);
                baseEEG8.etc.processStep = 'Rewarming';
                baseEEG8.etc.hoursBin=1;
                baseEEG8.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG8.srate)])>0))
                    pop_saveset(baseEEG8,'filename',char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG8.etc.processStep])));
                else
                    baseEEG8=[];
                end
            case 2
                baseEEG16=eeg_emptyset;
                baseEEG16 = pop_select(baseEEG48,'point',[8*60*60*baseEEG48.srate+1 8*2*60*60*baseEEG48.srate],'channel',1:baseEEG48.nbchan);
                [flatlines{i},baseEEG16] = findFlatlines(baseEEG16,baseDir,2,2);
                baseEEG16.etc.processStep = 'Rewarming';
                baseEEG16.etc.hoursBin=2;
                baseEEG16.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG16.srate)])>0))
                    pop_saveset(baseEEG16,char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG16.etc.processStep])));
                else
                    baseEEG16=[];
                end
            case 3
                baseEEG24=eeg_emptyset;
                baseEEG24 = pop_select(baseEEG48,'point',[8*2*60*60*baseEEG48.srate+1 8*3*60*60*baseEEG48.srate],'channel',1:baseEEG48.nbchan);
                [flatlines{i},baseEEG24] = findFlatlines(baseEEG24,baseDir,3,2);
                baseEEG24.etc.processStep = 'Rewarming';
                baseEEG24.etc.hoursBin=3;
                baseEEG24.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG24.srate)])>0))
                    pop_saveset(baseEEG24,char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG24.etc.processStep])));
                else
                    baseEEG24=[];
                end
            case 4
                baseEEG32=eeg_emptyset;
                baseEEG32 = pop_select(baseEEG48,'point',[8*3*60*60*baseEEG48.srate+1 8*4*60*60*baseEEG48.srate],'channel',1:baseEEG48.nbchan);
                [flatlines{i},baseEEG32] = findFlatlines(baseEEG32,baseDir,4,2);
                baseEEG32.etc.processStep = 'Rewarming';
                baseEEG32.etc.hoursBin=4;
                baseEEG32.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG32.srate)])>0))
                    pop_saveset(baseEEG32,char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG32.etc.processStep])));
                else
                    baseEEG32=[];
                end
            case 5
                baseEEG40=eeg_emptyset;
                baseEEG40 = pop_select(baseEEG48,'point',[8*4*60*60*baseEEG48.srate+1 8*5*60*60*baseEEG48.srate],'channel',1:baseEEG48.nbchan);
                [flatlines{i},baseEEG40] = findFlatlines(baseEEG40,baseDir,5,2);
                baseEEG40.etc.processStep = 'Rewarming';
                baseEEG40.etc.hoursBin=5;
                baseEEG40.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG40.srate)])>0))
                    pop_saveset(baseEEG40,char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG40.etc.processStep])));
                else
                    baseEEG40=[];
                end
            case 6
                baseEEG48 = pop_select(baseEEG48,'point',[8*5*60*60*baseEEG48.srate+1 8*6*60*60*baseEEG48.srate],'channel',1:baseEEG48.nbchan);
                [flatlines{i},baseEEG48] = findFlatlines(baseEEG48,baseDir,6,2);
                baseEEG48.etc.processStep = 'Rewarming';
                baseEEG48.etc.hoursBin=6;
                baseEEG48.etc.originalFiles= {files.name};
                if ~(size(flatlines{i},1)==1) || ~(sum(find([flatlines{i}]==1)>0) && sum(find([flatlines{i}==(8*60*60*baseEEG48.srate)])>0))
                    pop_saveset(baseEEG48,char(fullfile(baseDir,outputDir,['Bin_' num2str(i) '_' baseEEG48.etc.processStep])));
                else
                    baseEEG48=[];
                end
        end
    end
    
    function [baseEEG48]=initializeWindows(srate,nbchan,chanlocs)
         baseEEG48 = eeg_emptyset();
 
        baseEEG48.srate = srate;
        baseEEG48.nbchan = nbchan;
        baseEEG48.chanlocs = chanlocs;
        baseEEG48.data(1:baseEEG48.nbchan,1:(48*60*60)*baseEEG48.srate)=0;
        baseEEG48 = eeg_checkset(baseEEG48);

    end
end
