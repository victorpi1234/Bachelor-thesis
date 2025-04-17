function [ct, cst] = PTV_Shaped_Modulator(ct, cst, modulatorDepth, yLocation, cutout_start, HU, width_margingsYXZ, structureColor,  varargin)
% Creates a modulator shaped to PTV dimensions with specified depth
% Inputs:
%   ct, cst: CT and CST structures
%   modulatorDepth: Depth (y-dimension) of modulator (mm)
%   yLocation: Vertical offset from center
%   HU: Hounsfield units for modulator material
%   cutout_start: Offset from top of modulator where PTV cutout begins (mm)
%                 (0 = starts exactly at top surface)
% width_margingsYXZ : Margin of modulator width, compared to PTV (mm)

if ~exist('modulatorDepth','var') || isempty(modulatorDepth)
    modulatorDepth = 20; % Default depth in mm
end

if ~exist('HU','var') || isempty(HU)
    HU = 132;
end

if ~exist('cutout_start','var') || isempty(cutout_start)
    cutout_start = 0; % Default starts at top surface
end

if ~exist('width_margingsYXZ','var') || isempty(width_margingsYXZ)
    width_margingsYXZ = [5 5 5];
end

if ~exist('structureColor','var') || isempty(structureColor)
    structureColor = [0 1 1];
end


disp('Box_modulator function called!');

%% Find PTV to determine width
ixPTV = find(strcmp(cst(:,2), 'PTV'));
if isempty(ixPTV)
    error('PTV structure not found in CST');
end

% Get PTV voxels and calculate dimensions
PTV_voxels = cst{ixPTV,4}{1};
[ptv_y, ptv_x, ptv_z] = ind2sub(ct.cubeDim, PTV_voxels);

% Calculate PTV width with margin
x_margin = round(width_margingsYXZ(2)/ct.resolution.x);
z_margin = round(width_margingsYXZ(3)/ct.resolution.z); 
x_size = (max(ptv_x) - min(ptv_x)) + 2*x_margin;
z_size = (max(ptv_z) - min(ptv_z)) + 2*z_margin;
y_size = round(modulatorDepth/ct.resolution.y);

%% Create modulator structure
ixModulator = size(cst, 1) + 1; 
cst{ixModulator,1} = ixModulator-1;
cst{ixModulator,2} = 'modulator';
cst{ixModulator,3} = 'EXTERNAL';
cst{ixModulator,5} = struct('TissueClass',1, 'alphaX',0.1, 'betaX',0.05, 'Priority',2, 'Visible',1, 'visibleColor',[0 1 1], 'sfudOptimization',0);
ct.structureColor(ixModulator, :) = structureColor;

%% Create modulator volume
center = round([ct.cubeDim./2]);
Modulator_center = [center(1)-yLocation, center(2), center(3)];
halfSize = [y_size, x_size, z_size]/2;

% Create 3D mask
[Y,X,Z] = ndgrid(1:ct.cubeDim(1), 1:ct.cubeDim(2), 1:ct.cubeDim(3));
modulatorMask = (abs(X - Modulator_center(2)) <= halfSize(2)) & ...
               (abs(Y - Modulator_center(1)) <= halfSize(1)) & ...
               (abs(Z - Modulator_center(3)) <= halfSize(3));

%% Create PTV cutout starting from specified offset
% Calculate cutout parameters
cutout_start_vox = round(cutout_start/ct.resolution.y);
top_of_modulator = Modulator_center(1) - halfSize(1); % Top surface y-coordinate

% Move PTV voxels to start from cutout position
new_y = top_of_modulator + cutout_start_vox + (ptv_y - min(ptv_y));

% Keep only voxels that remain within bounds
valid = new_y >= 1 & new_y <= ct.cubeDim(1) & ...
        new_y <= (Modulator_center(1) + halfSize(1)); % Don't extend beyond bottom
moved_PTV = sub2ind(ct.cubeDim, new_y(valid), ptv_x(valid), ptv_z(valid));

% Cut out moved PTV from modulator
modulatorMask(moved_PTV) = false;

%% Finalize modulator
cst{ixModulator,4}{1} = find(modulatorMask);
ct.cubeHU{1}(cst{ixModulator,4}{1}) = HU;
end
