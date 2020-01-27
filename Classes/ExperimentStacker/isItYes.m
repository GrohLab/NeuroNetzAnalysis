function iOk = isItYes(str)
iOk = true;
if ~strcmpi(str,'y')
    fprintf(1,'Reply %s interpreted as ''no''\n', str);
    fprintf(1,'Quitting the object creation.\n');
    iOk = false;
end
end