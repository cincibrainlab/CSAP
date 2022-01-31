clear;
participant = 0; %PARTICIPANT NUMBER%
method='METHOD'; %PLACE METHOD (OToole or AB_Dynamic) HERE%
baseDir = 'INPUT DIR HERE'; %DIRECTORY OF PREPROCESSED INPUT FILES%
outputDir = fullfile(['OUTPUT DIR HERE' num2str(participant)]); %PLACE OUTPUT DIRECTORY IN BETWEEN QUOTES%
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
files = dir(fullfile(baseDir,'**/*.set'));
results={};
count=1;
absResults={};
absCount=1;
dbwpliResults = {};
dbwpliCount=1;
format long;

for m=1:length(files)
    EEG = pop_loadset(files(m).name,files(m).folder); 
    binName = char(regexp(EEG.setname,'Bin_\d*','match'));
    if contains(regexp(files(m).folder,'\w*.$','match'),{'followup','Non_Cooling'})
        processName = [char(regexp(EEG.setname,'(?<=\Bin_\d_)\w*(?=\_\d)','match')) '_' char(regexp(files(m).folder,'\d*.$','match'))];
    else
        processName = char(regexp(EEG.setname,'(?<=\Bin_\d_)\w*(?=\_\d)','match'));
    end

    %% 5 features 

    freq = 2:.5:64; 
    F = [2 3.5 4 7.5 8 10 10.5 12.5 13 30 30.5 55]; % 6-band
    [~, Locb] = ismember(F, freq); 
    Fidx = reshape(Locb, 2, []); nbband = size(Fidx,2); 

    srate = EEG.srate; 

    nepoch = EEG.trials; 

    nbchan = EEG.nbchan; 



    tic 

    absolute = NaN(nbchan, nepoch, nbband); ent = NaN(nbchan,nepoch);  

    flt = NaN(nbchan,nepoch); medianFreq = NaN(nbchan,nepoch);  

    hfd = NaN(nbchan,nepoch); 

    for k=1:EEG.nbchan 

        for i=1:EEG.trials 

            % 1st channel 

            epoch = EEG.data(k,:,i); % 16-sec 

            [~, ~, ~, P] = spectrogram(epoch, srate, round(srate/2), freq, srate); 

                     % window - 1s(srate), overlap - 50% (srate/2), F , Fs - srate 

        %     figure;subplot(3,1,[1 2]);imagesc(freq,[],P');ylabel('1-sec segments 50% overlapping') 

            P = mean(P,2); % averaging over windows 

        %     subplot(3,1,3);plot(freq,P);axis tight;xlabel('Frequency (Hz)') 

        %     saveas(gcf,['D:\Arc_Data\34\S01_FIGS\chan1_trial',num2str(i),'.png']);close 



            %%% median frequency (quantile) from absolute power spectrum based on 

            %%% 16-sec epoch 

            Pcum = cumsum(P); %figure;subplot(211);plot(freq,P);subplot(212);plot(freq,Pcum) 

            OverHalfIdx = find(Pcum>=Pcum(end)/2,1); 

            UnderHalfIdx = find(Pcum<=Pcum(end)/2,1,'last'); 

            MidFreq = (freq(OverHalfIdx)+freq(UnderHalfIdx))/2; 

            % fail case 1: power focused on 1st F (more than half of total) 

            if ~isempty(MidFreq) 

                medianFreq(k,i) = MidFreq;  

            else  

                medianFreq(k,i) = 1;  

            end 



            %%% band power (absolute) based on 16-sec epoch average 

            a = [mean(P(Fidx(1,1):Fidx(2,1))) mean(P(Fidx(1,2):Fidx(2,2))) mean(P(Fidx(1,3):Fidx(2,3))) mean(P(Fidx(1,4):Fidx(2,4))) mean(P(Fidx(1,5):Fidx(2,5))) mean(P(Fidx(1,6):Fidx(2,6)))]; 

            absolute(k,i,:) = a; 

            %%% >>> band ratios possible here 



            %%% relative band power (16-sec epoch) -> spectral entropy 

            p = a/sum(a); 

            ent(k,i) = -sum(p.*log(p))/log(6); 

            flt(k,i) = exp(1/6*sum(log(p)))/(1/6*sum(p)); 



            %%% Higuchi fractal dimension (16-sec epoch) 

            hfd(k,i) = Higuchi_FD(epoch, 10); % Kmax=10 

        end 

    end 

    toc % ~4.7 sec 

    pnt = Fidx';tic
    dat = gpuArray( double(EEG.data) );
    data = permute(dat, [2, 1, 3]); % timepnts*nchan*epochs
    datar = reshape(data, [], EEG.nbchan);
    Ftable = reshape(F, 2, [])';
    hil_complex = NaN(EEG.pnts, EEG.trials, EEG.nbchan, size(Ftable, 1));
    for i = 1:size(Ftable,1)    
        % fft gpu filter
        filt_order = round(9 * ( srate/Ftable(i, 1) ));
        d = designfilt('bandpassfir','FilterOrder', filt_order, ...
            'CutoffFrequency1',Ftable(i,1),'CutoffFrequency2',Ftable(i,2), ...
            'SampleRate', srate);
        B = d.Coefficients;
        filtered_data = fftfilt(gpuArray(B), datar);
        complex = zeros(EEG.pnts*nepoch, nbchan, 'gpuArray');
        for chani=1:nbchan
            complex(:, chani) = hilbert(filtered_data(:, chani));
        end     
        filtered_data = [];  
        hil_complex(:, :, :, i) = gather(reshape(complex, EEG.pnts, [], nbchan));    
        complex = [];
    end % 2048 x 150 x 6 x 6

    hcomplex = gpuArray( hil_complex ); % 1000x133x128x7
    dbwpli_mat = gpuArray(NaN(EEG.nbchan,EEG.nbchan, size(Ftable,1)));
    combos = combnk({EEG.chanlocs(:).labels}',2);
    for chancombo = 1 : length(combos)
        channel1 = combos{chancombo,1};
        channel2 = combos{chancombo,2};

        chan1 = find(strcmpi(channel1,{EEG.chanlocs.labels}));
        chan2 = find(strcmpi(channel2,{EEG.chanlocs.labels}));

        mass1 = squeeze(hcomplex(:,:,chan1,:)); % pnts*ntrial*band
        mass2 = squeeze(hcomplex(:,:,chan2,:));

        cdd_mass = mass1.*conj(mass2); % cross spectrum: (freq)pnts*ntrial*band
        dbwpli = (sum(imag(cdd_mass)).^2-sum(imag(cdd_mass).^2))./(sum(abs(imag(cdd_mass))).^2-sum(imag(cdd_mass).^2)); % 1*ntrials*bands
        dbwpli_mat(chan1,chan2,:) = squeeze(mean(squeeze(dbwpli))); % 6chan x 6chan x 6band
    end;toc % 1.115 sec

    %figure;for i=1:6;subplot(3,2,i);imagesc(squeeze(dbwpli_mat(:,:,i)));end

    %Features for accessing later 
    %(abberiviations based upon template,variables, and the full feature names)
    featureAbbreviations = {'medianFreq','ent','flt','hfd'}; 
    featureVars={medianFreq,ent,flt,hfd};
    features={'Median Frequency','Spectral entropy', 'Spectral flatness', 'Higuchi dimension'};
    featureDictionary = containers.Map(features,featureAbbreviations);
    columnNames = regexprep(regexprep(strip({EEG.chanlocs.labels}),' ',''),'-','_');

    for i=1:length(features) 
        lengthFeature = length(featureVars{i}); 
        for x=1:EEG.nbchan
            results{count,1} = EEG.chanlocs(x).labels;
            data = featureVars{i};
            results{count,2}=mean(data(x,:)');
            results{count,3}=participant;
            results{count,4}=string(processName);
            results{count,5}=string(binName);
            results{count,6} = string(featureDictionary(features{i}));
            count = count+1;
            data=[];
        end
    end
    
    bands={'delta', 'theta', 'alpha1', 'alpha2', 'beta','gamma'};
    for y=1:nbchan
        for j = 1:size(absolute,3)
            
            absResults{absCount,1} = EEG.chanlocs(y).labels;
            absResults{absCount,2}= bands{j};
            absResults{absCount,3}=mean(absolute(y,:,j)');
            absResults{absCount,4}=participant;
            absResults{absCount,5}=string(processName);
            absResults{absCount,6}=string(binName);
            absResults{absCount,7} = string('abs');
            absCount=absCount+1;
        end
    end
    
    for y=1:length(combos)
            
            dbwpliResults{dbwpliCount,1} = combos(y,1);
            dbwpliResults{dbwpliCount,2}= combos(y,2);
            dbwpliResults{dbwpliCount,3} = gather(dbwpli_mat(find(strcmpi(combos(y,1),{EEG.chanlocs.labels})),find(strcmpi(combos(y,2),{EEG.chanlocs.labels})),:));
            dbwpliResults{dbwpliCount,4}=participant;
            dbwpliResults{dbwpliCount,5}=string(processName);
            dbwpliResults{dbwpliCount,6}=string(binName);
            dbwpliResults{dbwpliCount,7} = string('dbwpli');
            dbwpliCount=dbwpliCount+1;
    end
    EEG=eeg_emptyset;
end

final=cell2table(results,'VariableNames',{'Channel','Value','ID','Process_Name','Bin_Name','Feature_Name'});
finalGrouped = groupsummary(final,["Channel" "ID" "Process_Name" "Bin_Name" "Feature_Name"],"mean");
writetable(finalGrouped,fullfile(outputDir,['Results' '.csv']));

finalAbs = cell2table(absResults,'VariableNames',{'Channel','Band','Value','ID','Process_Name','Bin_Name','Feature_Name'});
finalAbsGrouped = groupsummary(finalAbs,["Channel" "Band" "ID" "Process_Name" "Bin_Name" "Feature_Name"],"mean");
writetable(finalAbsGrouped,fullfile(outputDir,['Abs_Results' '.csv']));

finalDbwpli = cell2table(dbwpliResults,'VariableNames',{'Channel_1','Channel_2','Value','ID','Process_Name','Bin_Name','Feature_Name'});
writetable(finalDbwpli,fullfile(outputDir,['Dbwpli_Results' '.csv']));
temp=readtable(fullfile(outputDir,['Dbwpli_Results' '.csv']));
temp.Properties.VariableNames(3:length(bands)+2) = bands(:);
finalDbwpliGrouped = groupsummary(temp,["Channel_1" "Channel_2" "ID" "Process_Name" "Bin_Name" "Feature_Name"],"mean");
writetable(finalDbwpliGrouped,fullfile(outputDir,['Dbwpli_Results' '.csv']));

function [HFD] = Higuchi_FD(serie, Kmax) 
%{
Script for computing the Higuchi Fractal Dimension (HDF) of a signal.
INPUT:
    serie: is the temporal series that one wants to analyze by HDF. 
    It must be a row vector.
    Kmax: maximum number of sub-series composed from the original. To
    determine its values, we have followed the recommendation of Doyle et
    al at "Discriminating between elderly and young using a fractal 
    dimension analysis of centre of pressure". 
OUTPUT:
    HFD: the HFD of the temporal series.
PROJECT: Research Master in signal theory and bioengineering - University of Valladolid
DATE: 02/03/2014
AUTHOR: Jes�s Monge �lvarez
%}
%% Checking the ipunt parameters:
control = ~isempty(serie);
assert(control,'The user must introduce a series (first inpunt).');
control = ~isempty(Kmax);
assert(control,'The user must introduce the Kmax parameter (second inpunt).');
%% Processing:
% Composing of the sub-series:
N = length(serie); 
X = NaN(Kmax,Kmax,N);
for k = 1:Kmax
    for m = 1:k
        limit = floor((N-m)/k);
        j = 1;
        for i = m:k:(m + (limit*k))
            X(k,m,j) = serie(i);
            j = j + 1;
        end  
    end
end
% Computing the length of each sub-serie:
L = NaN(1, Kmax);
for k = 1:Kmax
    L_m = zeros(1,k);
    for m = 1:k
        R = (N - 1)/(floor((N - m)/k) * k);
        aux = squeeze(X(k,m,logical(~isnan(X(k,m,:))))); %We get the sub-serie without the NaNs.
        for i = 1:(length(aux) - 1)
            L_m(m) = L_m(m) + abs(aux(i+1) - aux(i));
        end
        L_m(m) = (L_m(m) * R)/k;
    end
    L(k) = sum(L_m)/k;
end
% Finally, we compute the HFD:
x = 1./(1:Kmax);
aux = polyfit(log(x),log(L),1);
HFD = aux(1); %We only want the slope, not the independent term. 
end