function ACEfit_Par = PrepData(ACEfit_Par)
%
% Updata 'ACEfit_Par' with the input data information
%
%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/34ea883
%          2018-09-21 09:59:07 +0100

% Start of "ACEfit" computations; shuffle into expected order.
ACEfit_Par    = Reorder(ACEfit_Par);

% Check Data
Vs            = CheckData(ACEfit_Par);

% Load Data
% Return an nElm x nSubj data matrix, where nElm = prod(Vs.Dim).
[Y,YM]        = LoadData(Vs);

% Process Data
ACEfit_Par    = ProcData(ACEfit_Par,Y,YM);

ACEfit_Par.Vs = Vs;

return
