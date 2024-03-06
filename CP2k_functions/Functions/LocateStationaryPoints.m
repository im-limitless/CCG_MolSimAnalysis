function [MinimaIndx] = LocateStationaryPoints(f)

MinimaIndx = islocalmin(f);

return

