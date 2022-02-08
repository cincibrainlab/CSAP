%-------------------------------------------------------------------------------
% remove_artefacts: simple procedure to remove artefacts
%
% Syntax: data=remove_artefacts(data,ch_labels,Fs,data_ref,ch_refs)
%
% Inputs: 
%     data      - EEG data, in bipolar montage; size: N_channels x N
%     ch_labels - cell of bipolar channel labels, 
%                 e.g. {'C3-O1','C4-O2', 'F3-C3', 'F4-C4'}
%     Fs        - sampling frequency (in Hz)
%     data_ref  - EEG data, in referential  montage; size: (N_channels+1) x N
%     ch_refs   - cell of referential channel labels, 
%                 e.g. {'C3','C4','F3','F4'}
%
% Outputs: 
%     data - EEG data after processing, in bipolar montage, size: N_channels x N
%
% Example:
%     Fs=256;
%     data_st=gen_test_EEGdata(2*60,Fs,1);
%     N=size(data_st.eeg_data_ref,2);
% 
%     % simulate artefacts:
%     % 1. F3 not properly attached:
%     if3=find(strcmp(data_st.ch_labels_ref,'F3'));
%     data_st.eeg_data_ref(if3,:)=randn(1,N).*10;
%
%     % 2. electrode coupling between C4 and Cz
%     ic4=find(strcmp(data_st.ch_labels_ref,'C4'));
%     icz=find(strcmp(data_st.ch_labels_ref,'Cz'));
%     data_st.eeg_data_ref(icz,:)=data_st.eeg_data_ref(ic4,:)+randn(1,N).*5;
%
%     % re-generate bipolar montage:
%     [data_st.eeg_data,data_st.ch_labels] = ...
%        set_bi_montage(data_st.eeg_data_ref,data_st.ch_labels_ref, ...
%                                 data_st.ch_labels_bi);
%
%     % remove channels:
%     eeg_art=remove_artefacts(data_st.eeg_data,data_st.ch_labels,data_st.Fs, ...
%                         data_st.eeg_data_ref,data_st.ch_labels_ref);


% John M. O' Toole, University College Cork
% Started: 05-04-2016
%
% last update: Time-stamp: <2019-09-05 13:59:26 (otoolej)>
%-------------------------------------------------------------------------------
function [data,timepoints]=remove_artefacts(data,ch_labels,Fs,params,data_ref,ch_refs)
if(nargin<3), error('requires 3 input arguments.'); end
if(nargin<4 || isempty(data_ref)), data_ref=[]; end
if(nargin<5 || isempty(ch_refs)), ch_refs=[]; end


%neural_parameters;

DBverbose=0;

[N_channels,N]=size(data);

%---------------------------------------------------------------------
% 0. check in referential mode first; is there problem with one 
%    channel (e.g. Cz)
%---------------------------------------------------------------------


%---------------------------------------------------------------------
% 1. look for electrode coupling:
%---------------------------------------------------------------------

ichannels=1:size(data,1);




% all other artefacts are on a channel-by-channel basis:
irem=[];
timepoints=struct([]);
for n=ichannels
    [data(n,:),timepoints]=art_per_channel(data(n,:),Fs,n,params,timepoints);
    
    irem=[irem find(isnan(data(n,:)))];
end

% remove artefacts across all channels
data(ichannels,unique(irem))=NaN;



function [x,timepoints]=art_per_channel(x,Fs,n,params,timepoints)
%---------------------------------------------------------------------
% remove artefacts on a per-channel basis
%---------------------------------------------------------------------
%neural_parameters;
DBverbose=1;

N=length(x);

%---------------------------------------------------------------------
% 1. electrode-checks (continuous row of zeros)
%---------------------------------------------------------------------
x_channel=x;
x_channel(x_channel~=0)=1;
irem=zeros(1,N);    
[lens,istart,iend]=len_cont_zeros(x_channel,0);
ielec=find(lens>=(params.ART_ELEC_CHECK*Fs));
timepoints(n).elec_cont_zeros_art = [istart(:), iend(:)];
if(~isempty(ielec) && DBverbose)
    fprintf(['electrode check at time(s): ' repmat('%g ',1,length(ielec))  ...
             ' \n'],istart(ielec)./Fs);
end
for m=ielec
    irun=[istart(m)-1:iend(m)+1];
    irun(irun<1)=1; irun(irun>=N)=N;
    irem(irun)=1;
    x(irun)=NaN;
end
timepoints(n).perc_elec_cont_zeros_art = 100*length(find(irem==1))/length(x);
if(any(irem==1) && DBverbose)
    fprintf('continuous row of zeros: %.2f%%\n', ...
            100*length(find(irem==1))/length(x));
end

x_nofilt=x;
[x_filt,inans]=filter_butterworth_withnans(x,Fs,50,0.5,[5 2]);


%---------------------------------------------------------------------
% 2. high-amplitude artefacts
%---------------------------------------------------------------------
art_coll=params.ART_TIME_COLLAR*Fs;
irem=zeros(1,N);    

x_hilbert=abs( hilbert(x_filt) );    

thres_upper=params.ART_HIGH_VOLT;
ihigh=find(x_hilbert>thres_upper);
timepoints(n).high_amp_art = ihigh';
if(~isempty(ihigh))
    for p=1:length(ihigh)
        irun=(ihigh(p)-art_coll):(ihigh(p)+art_coll);
        irun(irun<1)=1;  irun(irun>N)=N;               
        irem(irun)=1;
    end
end
x(irem==1)=NaN;
timepoints(n).perc_high_amp_art = 100*length(find(irem==1))/length(x);
if(any(irem==1) && DBverbose)
    fprintf('length of high-amplitude artefacts: %.2f%%\n', ...
            100*length(find(irem==1))/length(x));
end




%---------------------------------------------------------------------
% 3. continuous constant values (i.e. artefacts)
%---------------------------------------------------------------------
art_coll=params.ART_DIFF_TIME_COLLAR*Fs;
x_diff_all=zeros(1,N);
irem=zeros(1,N);    

x_diff_all=[diff(x) 0];        
x_diff=x_diff_all;
x_diff(x_diff~=0)=1;
[lens,istart,iend]=len_cont_zeros(x_diff,0);

% if exactly constant for longer than . then remove:
ielec=find(lens>=(params.ART_DIFF_MIN_TIME*Fs));
timepoints(n).const_values_art = [istart(:),iend(:)];

for m=ielec
    irun=[(istart(m)-art_coll):(iend(m)+art_coll)];
    irun(irun<1)=1; irun(irun>=N)=N;
    irem(irun)=1;
    x(irun)=NaN;
end

timepoints(n).perc_const_values_art = 100*length(find(irem==1))/length(x);
if(any(irem==1) && DBverbose)
    fprintf('continuous row of constant values: %.2f%%\n', ...
            100*length(find(irem==1))/length(x));
end


%---------------------------------------------------------------------
% 4. sudden jumps in amplitudes or constant values (i.e. artefacts)
%---------------------------------------------------------------------
art_coll=params.ART_DIFF_TIME_COLLAR*Fs;
irem=zeros(1,N);    
x_diff=x_diff_all;    
    
ihigh=find(abs(x_diff)>params.ART_DIFF_VOLT);
timepoints(n).sudden_jumps_art = ihigh';
if(~isempty(ihigh))
    for p=1:length(ihigh)
        irun=(ihigh(p)-art_coll):(ihigh(p)+art_coll);
        irun(irun<1)=1;  irun(irun>N)=N;               
        irem(irun)=1;
    end
end
xb=x;
x(irem==1)=NaN;


% before filtering, but should be eliminated anyway
x(inans)=NaN;
inans=find(isnan(x));
x_nofilt(inans)=NaN;
x=x_nofilt;

timepoints(n).perc_sudden_jump_art = 100*length(find(irem==1))/length(x);
if(any(irem==1) && DBverbose)
    fprintf('length of sudden-jump artefacts: %.2f%%\n', ...
            100*length(find(irem==1))/length(x));
    

end

