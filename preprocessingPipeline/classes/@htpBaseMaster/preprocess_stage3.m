%% Preprocess Stage 3
%
%%% Usage
%   
% obj = preprocess_stage3( obj )
%
%%% Parameters
%
% * INPUTS: obj
%
% * OUTPUTS: obj
%
% The optional input if the function is not self-invoked obj is the 
% htpPreprocessMaster object.  The output is the updated 
% htpPreprocessMaster object with information updated regarding the 
% status and details of stage 3 completion.
%
%%% Copyright and Contact Information
% Copyright (C) 2020 Cincinnati Children's (Pedapati Lab)
%
% This file is part of High Throughput Pipeline (HTP)
% 
% See https://bitbucket.org/eped1745/htp_stable/src/master/
% 
% Contact: ernest.pedapati@cchmc.org
%
% Revision 8/27 for better memory handling for large datasets


%% preprocess_stage3
% Performing principal component analysis to reduce the dimension and thus
% reduce amount of data and computation time.  PCA allows the data rank to
% be set and perform pre-reduction to hand off results to ICA computation.  
% ICA is performed to seperate various activity of sources and cortical areas
% through reduction of mutual information of channels.  
% The components produced by the PCA and ICA performance are then used for 
% preprocessing to locate the primary brain components and remove artifacts
% in later stages. 
%
% If the subject is on the incorrect stage, then that subject will be
% skipped in this stage of preprocessing.
%
% The user can specify which ICA algorithm to run for their computation of
% component analysis.  It is recommended to use the 'BINICA' option if the 
% user has the BINICA executable on their MATLAB path as this method is the 
% fastest computational route of those available in the pipeline.
function obj = preprocess_stage3( obj )

[mc, mm, mw] = obj.tools_log;

mm('Starting Stage 3: Final Cleaning...');

icadefs;  % from EEGLAB, must point to correct binica
mm(sprintf('ICA Processor: %s', ICABINARY));
mm(sprintf('Current Data Directory: %s', obj.htpcfg.basePath));

stage_last = 'preica';
stage_next = 'postica';

obj.sub = obj.loadSub( stage_last );

opt     = obj.formatOptions;

arrayfun(@( s ) s.setopt( opt ), obj.sub, 'uni', 0);

objStageStatus = obj.htpcfg.objStageStatus;
objStageStatus_completed = obj.htpcfg.objStageStatus_completed;

% obj.htpcfg.objStageStatus = find(obj.selectObjects(stage_last, obj.htpcfg.csvfile)); % current stage
% obj.htpcfg.objStageStatus_completed = find(obj.selectObjects(stage_next, obj.htpcfg.csvfile)); % completed


prev_files = 0; skip_files = 0; errorchk = 0;

sub = obj.sub;
     
savesub = @(x) obj.tool_saveSubject(x, stage_next);

%Switch for different cleaning methods
qEEG_process = 0;

for i = 1 : length(obj.sub)
    s=sub(i);

    if ismember(i, objStageStatus) % processes only files specified by the spreadsheet

        if strcmp('No', opt.always_recalc)
            try
                s.loadDataset( stage_next );
                s.proc_state = 'PostICA';
            catch
            end
        end

        s.loadDataset( 'preica' );
        
        if ~qEEG_process
            [abParams, s.EEG] = abDynamicParamInit(s.EEG);
            obj.msgout(sprintf('Parameters for AB cleaning:\nApproach: %s\nThreshold: %d microvolts\nFs: %d Hz\nWindow Size: %d s\n',abParams.Approach,...
                    abParams.Threshold,abParams.Fs,abParams.WindowSize),'step_complete');
            trials=size(s.EEG.data,3);
            s.EEG = abPrep(s.EEG,abParams,trials);
            s.EEG.data = reshape(s.EEG.data,s.EEG.nbchan,[],trials);
            
            obj.msgout(sprintf('AB'), 'step_complete');

            s.EEG = abPrep(s.EEG,abParams,trials);
            s.EEG.data = reshape(s.EEG.data,s.EEG.nbchan,[],trials);
            s.EEG = eeg_checkset(s.EEG);
        else
            s.autobadsegmentspreterm;
            obj.msgout(sprintf('OToole'), 'step_complete');
        end
        
        rejTrials = squeeze(any(sum(isnan(s.EEG.data))));
        if ~all(rejTrials)
            s.EEG = pop_rejepoch(s.EEG,rejTrials,0);
            s.epoch_badtrials = length(find(rejTrials));
            s.epoch_badid = ['[' num2str(find(rejTrials')) ']'];
            s.epoch_percent = 100-(s.epoch_badtrials / s.epoch_trials)*100;
            s.EEG = eeg_checkset(s.EEG);
            s.EEG.filename = s.filename.(stage_next);
            s.proc_state = 'PostICA';
            s.outputRow('postica');
            prev_files = prev_files + 1;


            % save cleaned dataset into the postica directory
            savesub(s);

            % unload data & decrease memory footprint
            s.unloadDataset;
        else
            if ismember(i, objStageStatus_completed)
                % unload data & decrease memory footprint
                s.unloadDataset;
                s.outputRow('postica');
                prev_files = prev_files + 1;
            else
                % unload data & decrease memory footprint
                s.unloadDataset;
                s.proc_state = 'ICA_Error';
                s.outputRow('error');
                skip_files = skip_files + 1;
            end
        end


    else  % if object is not at correct stage, push through only as saved

        if ismember(i, objStageStatus_completed)
            % unload data & decrease memory footprint
            s.unloadDataset;
            s.outputRow('postica');
            prev_files = prev_files + 1;
        else
            % unload data & decrease memory footprint
            s.unloadDataset;
            s.proc_state = 'ICA_Error';
            s.outputRow('error');
            skip_files = skip_files + 1;
        end
    end

    %obj.sub(i) = s;
    sub(i)=s;

end
    
obj.sub = sub;
sub=[];

% creates results spreadsheet
obj.msgout(sprintf('\nSummary'), 'step_complete');
obj.msgout(sprintf('\nFiles Processed: %d', length(objStageStatus)), 'step_msg');
obj.msgout(sprintf('Previously Completed Files: %d', prev_files), 'step_msg');
obj.msgout(sprintf('Skipped Files: %d', skip_files), 'step_msg');
obj.msgout(sprintf('Other Errors: %d', errorchk), 'step_msg');


obj.createResultsCsv(obj.sub, 'postica', 'Default');

end

% function [params,EEG] = abParamInit(EEG)
%     params.Approach = 'Total';
%     params.Threshold = 50;
%     params.Fs = EEG.srate;
%     params.WindowSize = 10; % unit in second
%     params.InData = double(EEG.data(reshape(find(~isnan(EEG.data(:,:,:))),size(EEG.data(:,:,:),1),[]))); % data matrix
%     %params.InData = double(reshape(EEG.data,EEG.nbchan,[]));
%     EEG.etc.abParams = struct();
%     EEG.etc.abParams.Approach = params.Approach;
%     EEG.etc.abParams.Threshold = params.Threshold;
%     EEG.etc.abParams.Fs = params.Fs;
%     EEG.etc.abParams.WindowSize = params.WindowSize;
%     EEG.etc.stage3Method = 'AB Non-Dynamic';
%     
% end

function [params,EEG] = abDynamicParamInit(EEG)
    prctile_threshold = 99.7;
    params.Approach = 'Window';
    params.Threshold = max(prctile(reshape(EEG.data,EEG.nbchan,[]),prctile_threshold,[2]));
    params.Fs = EEG.srate;
    params.WindowSize = 10; % unit in second
    params.InData = double(reshape(EEG.data,EEG.nbchan,[]));
    EEG.etc.abParams = struct();
    EEG.etc.abParams.Approach = params.Approach;
    EEG.etc.abParams.Prctile_Threshold = prctile_threshold;
    EEG.etc.abParams.Threshold = params.Threshold;
    EEG.etc.abParams.Fs = params.Fs;
    EEG.etc.abParams.WindowSize = params.WindowSize;
    EEG.etc.stage3Method = 'AB Dynamic';
end