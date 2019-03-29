MrToniDBPath = 'E:\Database\EphysData';
import2stack(MrToniDBPath)
[configStruct, expFile] = createConfigStruct(getParentDir(MrToniDBPath,1));
savePopConfigFile('ToniExample.gcf',configStruct)
[dPopStruct, cPopStruct, conditionStruct, configStruct] = getPopulationStack(getParentDir(MrToniDBPath,1));
PSTHStruct = getPopPSTH(dPopStruct,conditionStruct,configStruct);
plotPopPSTH(PSTHStruct,'BinSize',0.015,'PlotStyle','middle')
