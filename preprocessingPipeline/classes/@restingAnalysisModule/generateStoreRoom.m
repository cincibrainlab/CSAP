 function obj = generateStoreRoom( obj )

            EEG = obj.EEG; 
            dat = EEG.data;
            data = permute(dat, [2, 1, 3]); % pnts*nchan*trials(nseg)
            datap = permute(data, [1, 3, 2]); % pnts*trials(nseg)*nchan
            
            
            datapd = detrend3( datap );
            complex = zeros(size(datap));
            complex = fft( datap .* hann(length(datapd)));

            
            % absolute & normalized power
            EEG.abs_power = 2*abs(complex(1:EEG.pnts/2,:,:)).^2/(EEG.pnts*EEG.srate); % uV^2/Hz
            
            pnt = obj.pntsTable;n = size(pnt,1);
            intervals = []; 
            for k = 1:n
                intervals = [intervals, pnt(k,1):pnt(k,2)];
            end
            freq = obj.freqTable;
            EEG.rel_power = NaN*ones(size(EEG.abs_power));EEG.rel_power = EEG.rel_power(1:pnt(end,2),:,:);
            for i = 1:EEG.nbchan
                for j = 1:EEG.trials 


% relative power from here:
                    EEG.rel_power_sum(:, j, i) = EEG.abs_power(1:pnt(end,2), j, i)./ ...
                        sum(EEG.abs_power(intervals,j,i)); % align with rui_localpipeline
 
                    % denominator: sum over target frequency bands
                    
                end
            end
            


            obj.rest_abs_power = EEG.abs_power; % save as uV^2/Hz

            obj.rest_abs_hz = linspace(0, EEG.srate/2, EEG.pnts/2+1);
            obj.rest_rel_power = EEG.rel_power_sum;
            obj.rest_rel_hz = linspace(0, freq(end,2), pnt(end,2));

            
        end
        
       