function matrix_output=padded_cell2mat(cell_input,pad)
    % Convert a ragged cell array of arrays in a padded matrix
    % Pierre Morel 2015.
    
    
    if nargin<2
        pad=NaN;
    end

    %With for loop (faster than cellfun)
    
    lengths=cellfun(@numel,cell_input);
    maxLength=max(lengths);
    matrix_output=zeros(length(cell_input),maxLength)+pad;
    for k=1:length(cell_input)
        matrix_output(k,1:lengths(k))=cell_input{k};
    end
    
end