
function Goodfiles=matchfiles(thisdir,str,exclude)
Goodfiles={};
[goodfiles folders]=matchfiles_helper(thisdir,str,exclude);
Goodfiles=[Goodfiles;goodfiles];
while ~isempty(folders)
  thisfolder=folders{1};
  [goodfiles newfolders]=matchfiles_helper(thisfolder,str,exclude);
  Goodfiles=[Goodfiles;goodfiles];
  folders=[folders(2:end);newfolders];
end
    
