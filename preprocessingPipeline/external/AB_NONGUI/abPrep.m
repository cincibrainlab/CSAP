function [EEG_AB]=abPrep(EEG,abParam,trials)
%ABPREP Summary of this function goes here
%   Detailed explanation goes here
%     EEG = pop_eegfiltnew(EEG, .5, str2num(abParam.LowFilter), [], false, [], 0); % lowpass faster than bandpass
%     pop_eegplot(EEG,1,0,0)
    %pop_eegplot(EEG,1,0,0)
    EEG_AB = EEG;
    %indices = find(~isnan(EEG.data(:,:,:)));
     %zeroIndices = find(~isnan(abParam.InData));
%      EEG.data(zeroIndices) = 0;
     
    %abParam.InData = EEG.data(reshape(find(~isnan(EEG.data(:,:))),size(EEG.data(:,:),1),[])); % data matrix
    %abParam.InData(zeroIndices)=1;
    %abParam.InData = reshape(EEG.data(indices),[],1);
    tic;Parameters = Run_AB_rui(abParam);toc % quick
    %EEG_AB.data(indices) = reshape(Parameters.OutData,EEG.nbchan,[],trials);
    %EEG_AB.data(indices) = Parameters.OutData;
    EEG_AB.data = Parameters.OutData;
    EEG_AB.etc.abWarnings = Parameters.WarningsCount;
    %EEG_AB.data = reshape(EEG_AB.data,EEG.nbchan,[],trials);
    %pop_eegplot(EEG_AB,1,0,0)

    %% 
%     figure
%     for i=1:EEG.nbchan
%         subplot(2,3,i);h=histogram(Parameters.InData(i,:),-100:10:150);
%         h.Normalization = 'countdensity';xlabel('amplitude (muV)');ylabel('count');
%         title(EEG.chanlocs(i).labels)
%     end
%      suptitle(sprintf('Bin %s Step: %s with just non-NaN values',string(EEG_AB.etc.hoursBin), EEG_AB.etc.processStep))
%     subplot(2,3,1);h=histogram(Parameters.InData(1,:),-100:10:150);
%     h.Normalization = 'countdensity';xlabel('amplitude (muV)');ylabel('count');
%     title(EEG.chanlocs(2).labels)
end

