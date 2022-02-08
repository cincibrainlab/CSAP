function dh = detailHandler( obj, idx, operation, varargin )

if nargin > 3
    searchstr = varargin{1};
end

switch operation
    
    case 'init'
        
        strout = sprintf('Loading detail for obj %s', obj.sub(idx).subj_basename);
        obj.msgout(strout, 'msg_complete');
        
        dh = prepareDetails( obj.sub(idx) );

    case 'get' 
        
        dh = prepareDetails( obj.sub(idx) );
        
        dh = retrieveDetails( searchstr, dh.fn, obj.sub(idx) );

        
    otherwise
        obj.msgout('Invalid operation passed to detailHandler.', 'step_error');
end

    function dh = prepareDetails( s )
        
        dh.fn = fields( s );
        
        dh.basename =  s.subj_basename;
        dh.exclude = s.get_exclude_switch;
        dh.postcomp = checkifpostcomps( s );
        
    end

    function str = retrieveDetails( searchstr, fn, s )
        
        idx = find(strcmp(searchstr, fn));
        str = evalc('disp(s.(fn{idx}));');
        
    end

    function result = checkifpostcomps( sub )
        
        if strcmp('postcomps', sub.proc_state)
            
            result = true;
        else
            result = false;
        end
    end




end