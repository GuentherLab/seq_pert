function updatePermissions(filepath)
%
%   Updates permissions for a certain file or folder
%
%   INPUT:  filepath     Location of file or folder to be updated
%

if ismac
    permissioncall = ['chmod -R 770 ', filepath];
    [status, cmdout] = system(permissioncall);
    fprintf('\n=================\n\nFile permissions updated\n');
end
