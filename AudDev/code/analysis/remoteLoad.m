function variable = remoteLoad(remote, remoteFile)

% Function that loads a file from the scc

[~,~,ext] = fileparts(remoteFile);

if strcmp(ext, '.mat')
    
    if remote
        localFile = conn_cache('pull', remoteFile);
        variable = load(localFile);
    else
        variable = load(remoteFile);
    end
    
elseif strcmp(ext, '.csv')
    
    if remote
        localFile = conn_cache('pull', remoteFile);
        variable = readtable(localFile);
    else
        variable = readtable(remoteFile);
    end
end
