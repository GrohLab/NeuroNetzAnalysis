function [configStruct] = createConfigStruct(ExpDB, RecDB)
%CREATECOFIGSTRUCT returns a structure with the information required to
%proceed with a population analysis.
%   Prompting the user for specific inputs. For example, the cell type can
%   be red from the RecDB variable. However, the conditional variables or
%   signals need to be assumed and hard-coded in this function. 
if istable(RecDB) && istable(ExpDB)
    cellTypes = unique(RecDB.PhysioNucleus);
    []listdlg(
end

end

