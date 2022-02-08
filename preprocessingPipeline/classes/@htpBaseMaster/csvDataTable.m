function obj = csvDataTable( obj, stage, action )
% save and load

% selecting objects by CSV
stage_last = stage;
stage_next = 'NA';

switch action
    
    case 'load'
        
        load(obj.htpcfg.matfile, 'sub');
        arrayfun(@(sub) sub.setCsv( obj.htpcfg.csvfile ), sub, 'UniformOutput',false );
        arrayfun(@(sub) sub.setMat( obj.htpcfg.matfile ), sub, 'UniformOutput',false );
        arrayfun(@(sub) sub.updatePaths( obj.htpcfg.basePath), sub, 'UniformOutput',false );
        
        for i = 1 : length(sub)
           
            sub(i).htpcfg = obj.htpcfg;
            
        end
        
       
        obj.sub = arrayfun(@(x) x.outputRow(stage_last), sub);
        obj.createResultsCsv(obj.sub, stage_last, 'same');
        obj.dataTable = readtable(obj.htpcfg.csvfile);
        obj.sub = sub;
        
        
    case 'save'
        
        proc_stage_idx = find(strcmp('proc_state', obj.sub(1).log_subjHeader));
        new_proc_state_idx = find(strcmp('proc_state', obj.dataTable.Properties.VariableNames));        
        
        for i = 1 : length( obj.sub )
                       
            obj.sub(i).proc_state = char(obj.dataTable{i, new_proc_state_idx});
                    
        end
        
        arrayfun(@(x) x.outputRow(stage_last), obj.sub);
        try
        prompt = {'Create a suffix for this new dataset?'};
            dlgtitle = 'Input';
            dims = [1 35];
            definput = {'new_subset'};
            answer = char(inputdlg(prompt,dlgtitle,dims,definput));
        catch
            answer = 'new_subset';
        end
        
        obj.createResultsCsv(obj.sub, stage_last, answer );
        
        
        
        
        
end







end