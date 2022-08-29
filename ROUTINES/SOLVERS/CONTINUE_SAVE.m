function [filename] = CONTINUE_SAVE(U, J0, dUdlam, filename, backupname)
% Function for saving using CONTINUE.m call back feature

    if isfile(filename) || isfile([filename '.mat'])
        % File exists.
         
        % Load previous data
        tmp = load(filename);
    
        % Save backup in case of poor time out of slurm
        save(backupname, 'tmp');
        
        % Update data to save
        U = [tmp.U, U];
        dUdlam = cat(3, tmp.dUdlam, dUdlam);
        
    else
         % File does not exist.
    end
    
    save(filename, 'U', 'dUdlam');

end
