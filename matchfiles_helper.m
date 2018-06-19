
function [goodfiles folders]=matchfiles_helper(thisdir,str,exclude)
%thisdir=thisfolder
files=dir(thisdir);
%make list of files that match
goodfiles={};
folders={};
for i=1:numel(files)
    name=files(i).name;
    if numel(strfind(name,str))~=0,
        goodfiles=[goodfiles; [thisdir '\' name]];
    end
    if ~strcmp(name,'.') && ~strcmp(name,'..') && files(i).isdir
        excludeMatches=0;
        for ii=1:numel(exclude)
            excludeMatches=excludeMatches+strcmp(exclude{ii},name);
        end
        if excludeMatches==0
            folders=[folders; [thisdir '\' name]];
        end
    end
end