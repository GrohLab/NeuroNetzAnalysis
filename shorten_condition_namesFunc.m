function [condNames] = shorten_condition_names(condNames)
%Shorten Condition Names
%%Pattern for identifying frequency laser
CF="L"+digitsPattern(2)
%%Attachment for another condition modifier (e.g. DREADDS)
d="J60"
%%Renaming conditions into shortened format

for i=1:length(condNames)
    if contains(condNames(i),"Puff")
        if contains(condNames(i),"Control")
            if contains(condNames(i),d)
                condNames(i)="P"+"D"
            else
                condNames(i)="P"
            end
        elseif contains(condNames(i),"All")
            condNames(i)="PuffAll"
        end
    elseif contains(condNames(i),"Laser")
        if contains(condNames(i),"Control")
            if contains(condNames(i),d)
                condNames(i)="L"+"D"
            else
                condNames(i)="L"
            end
        elseif contains(condNames(i),"All")
            condNames(i)="LaserAll"
        end
    elseif contains(condNames(i),"Delay")
        if contains(condNames(i),CF)
            if sscanf(condNames(i),"Delay %f")>0.1
                if contains(condNames(i),d)
                    Fl=round(sscanf(condNames(i),"Delay %f"),1)*1000
                    condNames(i)=Fl+"F"+"D"
                else
                    Fl=round(sscanf(condNames(i),"Delay %f"),1)*1000
                    condNames(i)=Fl+"F"
                end
            elseif sscanf(condNames(i),"Delay %f")<0.1
                if contains(condNames(i),d)
                    Fb=round(sscanf(condNames(i),"Delay %f"),2)*1000
                    condNames(i)=Fb+"F"+"D"
                else
                    Fb=round(sscanf(condNames(i),"Delay %f"),2)*1000
                    condNames(i)=Fb+"F"
                end
            end
        else
            if sscanf(condNames(i),"Delay %f=")>0.1
                if contains(condNames(i),d)
                    Cl=round(sscanf(condNames(i),"Delay %f"),1)*1000
                    condNames(i)=Cl+"C"+"D"
                else
                    Cl=round(sscanf(condNames(i),"Delay %f"),1)*1000
                    condNames(i)=Cl+"C"
                end
            elseif sscanf(condNames(i),"Delay %f")<0.1
                if contains(condNames(i),d)
                    Cb=round(sscanf(condNames(i),"Delay %f"),2)*1000
                    condNames(i)=Cb+"C"+"D"
                else
                    Cb=round(sscanf(condNames(i),"Delay %f"),2)*1000
                    condNames(i)=Cb+"C"
                end
            end
        end
    end
end
end