function [EEGdata,counter] = AB_Correct(X,threshold)
% reference:
% Mourad 2007. A simple and fast algorithm for automatic suppression of
% high-amplitude artifacts in EEG data. 
%
% minimize the mean square error of Y-WX, where
% Y - target, W - smoothing mat (channel specific time weights)
% 
% J(W) = (Y-WX)' * (Y-WX)
% dJ(W)/dW = -2X' * (Y-WX) = 0 
% W = (XX')(YX')
% 
% X - input nbchan x time (default by outer function)
    counter=0;
    D = find(abs(X) >= threshold);
    Y = X;
    Y(D) = 0 ; % reference mat
    Rxx = X*X';
    Ryx = Y*X';
    
    lastwarn('');

    if sum(~isnan(Rxx))>=1
        W = Ryx * inv(Rxx); % smoothing mat
        EEGdata = W * X; % smoothed output
        [msgstr, msgid] = lastwarn;
        if strcmp(msgid,'MATLAB:illConditionedMatrix') || strcmp(msgid,'MATLAB:MATLAB:nearlySingularMatrix')
            counter=counter+1;
        end
    else
        EEGdata = X;
    end
end