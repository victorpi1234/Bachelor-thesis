function [ct] = makeCT(yDim, xDim, zDim, res)
disp('makeCT function called!');
ct.cubeDim      = [yDim xDim zDim]; % second cube dimension represents the x-coordinate

ct.resolution.x = res;
ct.resolution.y = res;
ct.resolution.z = res;
ct.numOfCtScen  = 1;

% create an ct image series with zeros - it will be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * -1000; % assign HU of Air
end