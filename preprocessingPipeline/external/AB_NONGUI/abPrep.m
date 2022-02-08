function [EEG_AB]=abPrep(EEG,abParam,trials)
%ABPREP Summary of this function goes here
%   Detailed explanation goes here

    EEG_AB = EEG;
   
     
    tic;Parameters = Run_AB_rui(abParam);toc % quick
    
    EEG_AB.data = Parameters.OutData;
    EEG_AB.etc.abWarnings = Parameters.WarningsCount;
    
end

