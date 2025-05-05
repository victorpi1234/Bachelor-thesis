function [ct, cst] = Box_modulator(ct, cst, boxSize, yLocation, removal_lenght, removal_spacing , tooth_height, HU , startpixel, varargin)
% boxSize in form [20, 130, 30] (how big the base for the modulator is)
% yLocation : where on the y-axis you want it to be (use 110 for example)



if ~exist('boxSize','var') || isempty(boxSize)
    boxSize = [20, 130, 30];
end
if ~exist('removal_lenght','var') || isempty(removal_lenght)
    removal_lenght = [2 2];
end
if ~exist('removal_spacing','var') || isempty(removal_spacing)
    removal_spacing = 2*removal_lenght;
end
if ~exist('tooth_height','var') || isempty(tooth_height)
    tooth_height = boxSize(1)/2;
end
if ~exist('HU','var') || isempty(HU)
    HU = 132;
end
if ~exist('startpixel','var') || isempty(startpixel)
    startpixel = removal_lenght;
end


if numel(removal_spacing) ~= 2
    error('removal_spacing must be a 2-element vector [x z]');
end
if numel(removal_lenght) ~= 2
    error('removal_lenght must be a 2-element vector [x z]');
end
% Warning
if boxSize(1) < tooth_height
    warning('Box height (%d) is smaller than tooth height (%d). Adjusting tooth height.', boxSize(1), tooth_height);
    tooth_height = boxSize(1) / 2;
end


space_between_removal_x = removal_spacing(1);
space_between_removal_z = removal_spacing(2);
lenght_of_removal_x = removal_lenght(1);
lenght_of_removal_z = removal_lenght(2);
startpixel_x = startpixel(1);
startpixel_z = startpixel(2);

% if nargin < 8 || isempty(HU)
%     HU = 300;
% end
% 

% p = inputParser;
% addOptional(p, 'HU', 300);
% addOptional(p, 'startpixel_x', 0);
% addOptional(p, 'startpixel_z', 0);
% addOptional(p, 'space_between_removal_x', 3);
% addOptional(p, 'space_between_removal_z', 3);
% addOptional(p, 'lenght_of_removal_x', 1);
% addOptional(p, 'lenght_of_removal_z', 1);
% addOptional(p, 'tooth_height', boxSize(1)/2);
% 
% parse(p, varargin{:});
% 
% HU = p.Results.HU;
% startpixel_x = p.Results.startpixel_x;
% startpixel_z = p.Results.startpixel_z;
% space_between_removal_x = p.Results.space_between_removal_x;
% space_between_removal_z = p.Results.space_between_removal_z;
% lenght_of_removal_x = p.Results.lenght_of_removal_x;
% lenght_of_removal_z = p.Results.lenght_of_removal_z;
% tooth_height = round(p.Results.tooth_height);


disp('Box_modulator function called!');
%% Error messages

if any(mod(boxSize, 2) ~= 0)
    error('Box_modulator:InvalidInput', 'boxSize must be even in all dimensions. Received [%s].', num2str(boxSize));
end

%%

% ixModulator = 4;

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

%%
% Extract PTV voxel indices and center
% Define the box size
%cst{ixModulator,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(5*10^3,11.6));
%boxSize = [20, 130, 30]; % Adjust as needed
center = round([ct.cubeDim./2]);
Modulator_center =[(center(1)- yLocation) center(2) center(3)];
% Create a modulator box around the PTV center
modulatorHelper = zeros(ct.cubeDim);
halfSize = boxSize / 2;



for x = 1:ct.cubeDim(2)
    for y = 1:ct.cubeDim(1)
        for z = 1:ct.cubeDim(3)
            if abs(x - Modulator_center(2)) <= halfSize(2) && ...
                abs(y - Modulator_center(1)) <= halfSize(1) && ...
                abs(z - Modulator_center(3)) <= halfSize(3)
                modulatorHelper(y,x,z) = 1;
            end
        end
    end
end

% Remove alternating x-slices starting from bottom up to tooth_height
for i = startpixel_x:space_between_removal_x:2*halfSize(2)  % Locate every second slice in x-direction
    x = Modulator_center(2) - halfSize(2) + i;  % x position from beginning to end of modulator
    if x > 0 && x <= ct.cubeDim(2) 
        % Start from bottom (Modulator_center(1) + halfSize(1)) and go up to tooth_height
        for y = (Modulator_center(1) + halfSize(1)):-1:(Modulator_center(1) + halfSize(1) - tooth_height)
            for z = 1:ct.cubeDim(3)
                if abs(z - Modulator_center(3)) <= halfSize(3)
                    %x_range = max(1,
                    %floor(x-lenght_of_removal_x/2)):min(ct.cubeDim(2),
                    %ceil(x+lenght_of_removal_x/2)); % this inputs bigger
                    %ranges for very small values
                    half = floor(lenght_of_removal_x / 2);
                    if mod(lenght_of_removal_x, 2) == 0  % even
                        x_range = max(1, x - half + 1):min(ct.cubeDim(2), x + half);
                    else  % odd
                        x_range = max(1, x - half):min(ct.cubeDim(2), x + half);
                    end
                    modulatorHelper(y, x_range, z) = 0; % Remove part of the modulator
                end
            end
        end
    end
end

% Remove alternating z-slices starting from bottom up to tooth_height
for i = startpixel_z:space_between_removal_z:2*halfSize(3)  % Locate every second slice in z-direction
    z = Modulator_center(3) - halfSize(3) + i;  % z position from beginning to end of modulator
    if z > 0 && z <= ct.cubeDim(3) 
        % Start from bottom (Modulator_center(1) + halfSize(1)) and go up to tooth_height
        for y = (Modulator_center(1) + halfSize(1)):-1:(Modulator_center(1) + halfSize(1) - tooth_height)
            for x = 1:ct.cubeDim(2)
                if abs(x - Modulator_center(2)) <= halfSize(2)
                    z_range = max(1, floor(z-lenght_of_removal_z/2)):min(ct.cubeDim(3),ceil(z+lenght_of_removal_z/2));
                    modulatorHelper(y, x, z_range) = 0; % Remove part of the modulator
                end
            end
        end
    end
end

% Store voxel indices for the modulator in `cst`
cst{ixModulator,4}{1} = find(modulatorHelper);
Ind_modulator = cst{ixModulator,4}{1};
%%
% Assign HU values in the CT image
ct.cubeHU{1}(cst{ixModulator,4}{1}) = HU; % Adjust HU value if needed
ct.cubeHU{1}(Ind_modulator) = HU; 
