function pln = specific_pln_maker(ct, cst, doseGridRes)
disp("Pln made!")

pln.radiationMode = 'protons';
pln.machine               = 'Generic';
pln.numOfFractions        = 1;
pln.propStf.gantryAngles  = [0];
pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% pln.propStf.longitudinalSpotSpacing = 8;
% pln.propStf.bixelWidth = 10;
% doseGridRes = 1;
pln.propDoseCalc.doseGrid.resolution = struct('x',doseGridRes,'y',doseGridRes,'z',doseGridRes); %[mm]

pln.bioModel = matRad_bioModel(pln.radiationMode,'none');
% pln.propOpt.bioOptimization = 'none';
pln.propDoseCalc.calcLET = 0;
pln.propDoseCalc.engine = 'HongPB';

end 