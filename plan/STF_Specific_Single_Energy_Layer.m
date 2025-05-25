function [stf] = STF_Specific_Single_Energy_Layer(ct, cst, pln)
disp('Result function called!');
%% Generate Beam Geometry STF
% pln.propStf.addMargin    = false; %to make smaller stf, les bixel
% create multi energy stf
mb_stf = matRad_generateStf(ct,cst,pln);
stf = mb_stf;
%%
% stf = mb_stf;
% stf_no_mod = mb_stf_no_mod;
%%
% create single bixel stf
sb_stf = matRad_generateSingleBixelStf(ct,cst,pln);
% sb_stf_no_mod = matRad_generateSingleBixelStf(ct_no_mod,cst_no_mod,pln);
%%
% adapt stf to have one energy layer/bixel per ray
% run through all rays, keep only ray in the (x, 0, 0) line
% to make a one front beam
% can be altered to have another line
list = []; % empty
for j = 1:stf.numOfRays
    % if j==111
    if (stf.ray(j).rayPos_bev(2)~=0 || stf.ray(j).rayPos_bev(3)~=0)
% tmp.ray(j) = [];
  list = [list; j]; % save all rows that have (x, 0, 0) coordinates
    end
    % end
end

f = flip(list); % flip to start removing from the bottom
stf.ray(f) = []; % delete all rows that have (x, 0, 0) coordinates

% remove fields because they interfeer with code
stf.ray = rmfield(stf.ray, 'numParticlesPerMU');
stf.ray = rmfield(stf.ray, 'minMU');
stf.ray = rmfield(stf.ray, 'maxMU');

% make energy, rangeShifter and focusIx the same as in sb_stf
for i = 1:size(stf.ray,2)
    %stf.ray(i).energy = sb_stf.ray(1).energy;
    %stf.ray(i).energy =1.184896818197490e+02;
    % stf.ray(i).energy =1.104896818197490e+02;
    stf.ray(i).energy =1.070896818197490e+02;
    stf.ray(i).rangeShifter = sb_stf.ray(1).rangeShifter;
    stf.ray(i).focusIx = sb_stf.ray(1).focusIx;
end

% change remainning parameters
stf.numOfRays  = size(stf.ray,2);
stf.numOfBixelsPerRay = ones(stf.numOfRays,1)';
stf.totalNumOfBixels = sum(stf.numOfBixelsPerRay(:));

%figure, matRad_plotSliceWrapper(gca,ct_comb,cst_comb,1,resultGUI.physicalDose - resultGUI_comb.physicalDose ,3,slice);
% slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
% slice = slice(3);

end