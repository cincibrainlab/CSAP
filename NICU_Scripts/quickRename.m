function quickRename(baseDir,varargin)
    if ~isempty(varargin)
        containers.Map({varargin{1}},{repelem('',length(varargin{1}))});
    end
    files = dir(fullfile(baseDir, '*.EDF'));
    metadata=cell(length(files),1);
    startInfo = strings(length(files),2);
    for i=1:length(files)
        [metadata{i},~] = import_edf(fullfile(baseDir,files(i).name),0);
        startInfo(i,1) = regexprep(datestr(regexp(metadata{i}.recordID,'\d+\-\w*\-\d+','match'),23),'/','-');
        startInfo(i,2) = [metadata{i}.startTime];
        startInfo(i,2) = regexprep(metadata{i}.startTime,'\.','-');
        newName=regexprep(regexp(regexprep(strip(metadata{i}.patientID),' ',''),'(?<=\d{2}\-\w*\-\d{4}).*$','match'),'\,','-');
        if isempty(newName)
            newName=regexprep(regexp(regexprep(strip(metadata{i}.patientID),' ',''),'(?<=\XXX).*$','match'),'\,','-');
        end
        newName = regexprep(newName,'\"','');
        newName = join([join([newName,startInfo(i,1),startInfo(i,2)],'_') 'EDF'],'.');
        if ~exist(fullfile(baseDir,newName),'file')
            movefile(fullfile(baseDir,files(i).name),fullfile(baseDir,newName));    %change the filenames
        else
            continue
        end
    end
end