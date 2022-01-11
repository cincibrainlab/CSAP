function [finalFlatlinePoints,EEG] = findFlatlines(EEG,baseDir,binIndex,hourSegment)
%CLEAN Summary of this function goes here
    %   Detailed explanation goes here
    max_flatline_duration = 60;
    numberPointsToBeCut = 0;
    if ~exist('max_flatline_duration','var') || isempty(max_flatline_duration) max_flatline_duration = 5; end
    if ~exist('max_allowed_jitter','var') || isempty(max_allowed_jitter) max_allowed_jitter = 20; end
    flatlineNum = zeros(EEG.nbchan,1);
    % flag segments
    for c=1:EEG.nbchan
        zero_intervals = reshape(find(diff([false abs(diff(EEG.data(c,:)))<(max_allowed_jitter*eps) false])),2,[])';            
        temp  = find(zero_intervals(:,2) - zero_intervals(:,1) > max_flatline_duration*EEG.srate);
        temp_length = length(temp);
        flatlines{c} =  [zero_intervals(temp,1) zero_intervals(temp,2);];
        flatlineNum(c) = temp_length;
    end
    if(~any(flatlineNum(1)==flatlineNum))
        error('There is an issue where one of the channels is not dead at the same time as the others\n');
    else
        finalFlatlinePoints = flatlines{1};
    end
    if ~isempty(finalFlatlinePoints) && (size(finalFlatlinePoints,1)>1 || ~(sum(find([finalFlatlinePoints(1,1)]==1)>0) && sum(find([finalFlatlinePoints(1,2)==(8*60*60*EEG.srate)])>0)))
        fig=figure('Name',['Bin_' num2str(binIndex) '_Cuts'],'units','normalized','outerposition',[0 0 1 1],'visible','off');
        %plot(EEG.times/1000/60/60,mean(EEG.data./EEG.srate./60./60,1),'black')
        plot(EEG.times/1000/60/60,mean(EEG.data,1),'black')
        title('Cuts to be Made to Pre-Filtered EEG Data (Averaged Across Channel)')
        xlabel('Time (Hours)')
        ylabel('Amplitude (µV)')
        %limitValue = max(max(mean(EEG.data./EEG.srate./60./60)))+.001;
        limitValue = max(max(abs(mean(EEG.data,1))));
        ylim([-limitValue,limitValue])
        %ylim([-limitValue,limitValue])
        y_limits = ylim;
        for i=1:size(finalFlatlinePoints,1)
            vertices = [EEG.times(finalFlatlinePoints(i,1))/1000/60/60 y_limits(2); EEG.times(finalFlatlinePoints(i,1))/1000/60/60 y_limits(1); EEG.times(finalFlatlinePoints(i,2))/1000/60/60 y_limits(2); EEG.times(finalFlatlinePoints(i,2))/1000/60/60 y_limits(1);];
            faces = [3 4 2 1;];
            patch('Faces', faces, 'Vertices', vertices, 'FaceColor', 'red','FaceAlpha', 0.1);
            EEG.data(:,finalFlatlinePoints(i,1):finalFlatlinePoints(i,2)) = NaN;
            numberPointsToBeCut = numberPointsToBeCut + numel(finalFlatlinePoints(i,1):finalFlatlinePoints(i,2));
        end
        EEG.etc.percentageFlatline = double((numberPointsToBeCut/size(EEG.data,2))*100);
        %EEG.etc.hourSegmentation = hourSegment;
        EEG.etc.flatlinePointRanges = {finalFlatlinePoints};
        %EEG.etc.existingDataIndices = find(~isnan(EEG.data(:,:)));
        if ~exist(fullfile(baseDir,'Cut_Images'),'dir')
            mkdir(fullfile(baseDir,'Cut_Images'));
        end
        saveas(fig,fullfile(baseDir,'Cut_Images',['Bin_' num2str(binIndex) '_Cuts.jpg']),'jpeg');
        close(['Bin_' num2str(binIndex) '_Cuts'])
    end
    EEG.etc.hourSegmentation = hourSegment;
end