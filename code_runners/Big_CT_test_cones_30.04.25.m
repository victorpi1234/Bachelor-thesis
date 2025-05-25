%%
clear all
close all
clc
%% set matRad configuration
addpath(genpath('C:\Users\victo\OneDrive\Работен плот\dskalo1\Heidelberg Physik\bachelor_thesis\dkfz\matRad-dev'));
matRad_rc;
%%
addpath(genpath('C:\Users\victo\OneDrive\Работен плот\dskalo1\Heidelberg Physik\bachelor_thesis\dkfz\MymatRad\Project'));
%%
ct_no_mod = makeCT(400,200,100, 1);
[ct_no_mod, cst_no_mod] = makeCST_water(ct_no_mod, 30, 20, [-50 0 0], [],[],[],-10);
% Check for resume state
% resume = false;
% if isfile('resume_state.mat')
%     load('resume_state.mat');  % Loads `currentState`
%     resume = true;
% end

%% stf


clearvars ct cst
[ct,cst] = new_PTV_shaped_coned_inverse_good(ct_no_mod,cst_no_mod,6, 25,3,200,6,6,132, [] ,[1 1 1]);
%%
pln = specific_pln_maker(ct, cst, 0.5);
%%
stf = calc_STF_Single_Energy_Layer_Beam(ct, cst, pln);


% Calculation without optimization
w = ones(stf.totalNumOfBixels,1);
%%
resultGUI = matRad_calcDoseForward(ct,cst,stf,pln,w);
resultGUI_no_mod = matRad_calcDoseForward(ct_no_mod,cst_no_mod, stf, pln, w);
%% Title and analysis
title_text = "namebob";
data_analysis_one_Dose(ct, cst, resultGUI, pln, title_text, 50, "on");

%% Bragg peak 9 pictures
ixProfileX = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
ixProfileX = ixProfileX(1)-8:ixProfileX(1)+1;
ixProfileY = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
for ixProfileY = ixProfileY(2):ixProfileY(2)+1
    slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
    slice = slice(3);
    profile = resultGUI.physicalDose(:,ixProfileY,slice);
    % Store both profiles for comparison
    
    figure,plot(profile,'LineWidth',2),grid on,hold on,
    plot(profile,'LineWidth',2),legend({'original profile'}),
    xlabel('mm'),ylabel('Gy(RBE)'),title('profile plot')
end

%% Save png
folderpng = 'C:\Users\victo\OneDrive\Работен плот\dskalo1\Heidelberg Physik\bachelor_thesis\dkfz\DKFZ_BachelorThesis\presentations\pictures\';
saveas(gcf, [folderpng,'Dose_in_y_direction_1.png'])


%%
matRadGUI

% To fix: boxSize_y doesn't work with odd numbers
% Make also, that the names of every file, meaning  title_text = sprintf('Cone-width=3cm Cone-spacing=5cm \n s-water-to-phantom = 35cm HU = 132');