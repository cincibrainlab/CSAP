
function AB_Parameters = Run_AB(Parameters)
%function AB_Parameters = Run_AB(Parameters)

AB_Parameters = Parameters


%[eeg] = read_eep_cnt(Parameters.InputFile,Parameters.FirstSamp,Parameters.LastSamp);  %reading the EEG data





%Running AB

fprintf('\n Running AB .....');


AB_alg_window_GUI();  %Running AB

fprintf('Done. \n')

%Writing the data to the output file


%=======================================================================
    function   AB_alg_window_GUI()
        %function  EEGdata = AB_alg_window_GUI(oldX,threshold,Approach);

        threshold = AB_Parameters.Threshold;
        [N,T] = size(AB_Parameters.InData);

        if AB_Parameters.Approach == 'Total '
            LWind = T;
        
        elseif AB_Parameters.Approach == 'Window'
            %LWind = min([T ,100*N]);
            LWind = min([T ,10*AB_Parameters.Fs]);
        else
            error('Unknown approach parameter');
        end


        %======================================================================
        % removing the mean value of the input data
        %======================================================================
        for i = 1:N
            meanValue = mean(AB_Parameters.InData(i,:));
            AB_Parameters.InData(i,:) = AB_Parameters.InData(i,:) - meanValue ;
        end
        %======================================================================
        % Removing the artifacts
        %======================================================================
        if LWind == T
            AB_Parameters.OutData =  AB_Parameters.InData;
            D = find(abs(AB_Parameters.OutData) >= threshold);
            AB_Parameters.OutData(D) = 0;
            Rxx = AB_Parameters.InData*AB_Parameters.InData';
            Ryx = AB_Parameters.OutData*AB_Parameters.InData';
            W = Ryx * inv(Rxx);
            AB_Parameters.OutData = W * AB_Parameters.InData;
        else
            AB_Parameters.OutData = zeros(N,T);
            Iin = 1; Ifin = LWind + 1; r = 0.02;
            while 1
                ind = 1;LWind2 = LWind;
                if Iin == 1; ind = 0; end
                if Ifin > T; ind = 2; Ifin = T; LWind2 = Ifin - Iin; end
                W = mywindow(LWind,r,ind,LWind2);
                W = W(:)';
                X = AB_Parameters.InData(:,Iin:Ifin -1);
                %     [size(X) size(W)]
                for i = 1:N

                    X(i,:) = W.*X(i,:);
                end
                %     [size(EEGdata(:,Iin:Ifin-1)) size(X)]
                AB_Parameters.OutData(:,Iin:Ifin-1) = AB_Parameters.OutData(:,Iin:Ifin-1) + AB_correct(X,threshold);
                if Ifin >= T; break; end
                Iin = Ifin - fix(r*LWind/4) + 1;
                Ifin = Iin + LWind;
                %     [T Ifin]
            end
        end
    end
end



%#####################################################################
function EEGdata = AB_correct(X,threshold);
[N,M]=  size(X);
if M > N
    D = find(abs(X) >= threshold);
    Y = X;
    Y(D) = 0 ;
    Rxx = X*X';
    Ryx = Y*X';
    W = Ryx * inv(Rxx);
    EEGdata = W * X;
    
else
    EEGdata = zeros(N,M);
end

end
