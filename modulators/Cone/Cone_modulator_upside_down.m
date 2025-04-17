function [ct, cst] = Cone_modulator_upside_down(ct, cst, coneBaseRadius, coneHeight, yLocation, xLocation, baseThickness, boxSize, coneSpacingX_Z, HU, varargin)
% boxSize in form [20, 130, 30] (how big the base for the modulator is)
% yLocation : where on the y-axis you want it to be

if ~exist('coneSpacingX_Z','var') || isempty(coneSpacingX_Z)
    coneSpacingX_Z = [2*coneBaseRadius 2*coneBaseRadius];
end
if ~exist('boxSize','var') || isempty(boxSize)
    boxSize = [20, 130, 30];
end
if ~exist('xLocation','var') || isempty(xLocation)
    xLocation = 0;
end
if ~exist('HU','var') || isempty(HU)
    HU = 132;
end

% if nargin < 3 || isempty(coneBaseRadius)
%     coneBaseRadius = 5;
% end
% if nargin < 4 || isempty(coneHeight)
%     coneHeight = 30;
% end
% % if nargin < 5 || isempty(boxSize)
% %     boxSize = [20, 130, 30];
% % end
% if nargin < 5 || isempty(yLocation)
%     yLocation = 110;
% end

% if nargin < 8 || isempty(HU)
%     HU = 300;
% end
p = inputParser;
% addOptional(p, 'HU', 300);
% addOptional(p, 'boxSize', [20, 130, 30]);
%addOptional(p, 'coneHeight', 30);
parse(p, varargin{:});
% HU = p.Results.HU;
% boxSize = p.Results.boxSize;
disp('Cone_modulator function called');


coneSpacingX = coneSpacingX_Z(1); 
coneSpacingZ = coneSpacingX_Z(2);

%%

ixModulator = size(cst, 1) + 1; 

% define general VOI properties
cst{ixModulator,1} = ixModulator-1;
cst{ixModulator,2} = 'modulator';
cst{ixModulator,3} = 'EXTERNAL';

% define optimization parameter for both VOIs
cst{ixModulator,5}.TissueClass  = 1;
cst{ixModulator,5}.alphaX       = 0.1000;
cst{ixModulator,5}.betaX        = 0.0500;
cst{ixModulator,5}.Priority     = 2;
cst{ixModulator,5}.Visible      = 1;
cst{ixModulator,5}.visibleColor = [0 1 1];
cst{ixModulator,5}.sfudOptimization=0;
ct.structureColor(ixModulator, :) = [0 1 1];

%% Define cone parameters
% coneSpacingX = 2*coneBaseRadius; % Spacing between cones along the x-axis
% coneSpacingZ = 2*coneBaseRadius; % Spacing between cones along the z-axis
center = round([ct.cubeDim./2]);

% Define the box size for the modulator
%boxSize = [20, 130, 30]; % Adjust as needed
Modulator_center = [(center(1) - yLocation), center(2) - xLocation, center(3)]; % Center of the modulator (center(1)-110) for original dimension

% Create a modulator box around the PTV center
modulatorHelper = zeros(ct.cubeDim);
halfSize = boxSize / 2;
%%
% Create Base (n pixels thick, on top of cones)
baseTopY = Modulator_center(1) + baseThickness - 1;  % Top of the base

% Ensure base does not exceed cube dimensions
if baseTopY > ct.cubeDim(1)
    baseTopY = ct.cubeDim(1);  % Clamp to max Y-dimension
end

% Fill the base (n pixels thick)
for x = 1:ct.cubeDim(2)
    for y = Modulator_center(1):baseTopY  % Only fill 'n' pixels upward
        for z = 1:ct.cubeDim(3)
            if abs(x - Modulator_center(2)) <= halfSize(2) && ...
               abs(z - Modulator_center(3)) <= halfSize(3)
                modulatorHelper(y,x,z) = 1;  % Fill the base
            end
        end
    end
end

% Loop to create conical teeth (pointing downward)
for x = Modulator_center(2) - halfSize(2):coneSpacingX:Modulator_center(2) + halfSize(2)
    for z = Modulator_center(3) - halfSize(3):coneSpacingZ:Modulator_center(3) + halfSize(3)
        % Create a cone at (x, y, z)
        for y = Modulator_center(1)-1:-1:Modulator_center(1) - coneHeight  % Grow downward
            % Calculate the current radius of the cone (tapered)
            currentRadius = coneBaseRadius * (1 - (Modulator_center(1) - y) / coneHeight);
            for dx = -ceil(currentRadius):ceil(currentRadius)
                for dz = -ceil(currentRadius):ceil(currentRadius)
                    % Check if the voxel is within the cone radius
                    if dx^2 + dz^2 <= currentRadius^2
                        % Calculate voxel indices
                        voxelX = x + dx;
                        voxelY = y;
                        voxelZ = z + dz;
                        % Ensure voxel indices are within bounds
                        if voxelX >= 1 && voxelX <= ct.cubeDim(2) && ...
                           voxelY >= 1 && voxelY <= ct.cubeDim(1) && ...
                           voxelZ >= 1 && voxelZ <= ct.cubeDim(3)
                            modulatorHelper(voxelY, voxelX, voxelZ) = 1;
                        end
                    end
                end
            end
        end
    end
end

cst{ixModulator,4}{1} = find(modulatorHelper);
Ind_modulator = cst{ixModulator,4}{1};

ct.cubeHU{1}(cst{ixModulator,4}{1}) = HU; % Adjust HU value if needed
ct.cubeHU{1}(Ind_modulator) = HU; 

end