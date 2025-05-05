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
ct_no_mod = makeCT(200,200,100, 1);

% Check for resume state
% resume = false;
% if isfile('resume_state.mat')
%     load('resume_state.mat');  % Loads `currentState`
%     resume = true;
% end


for PTV_center = [-50]
    [ct_no_mod, cst_no_mod] = makeCST_water(ct_no_mod, 30, 20, [PTV_center 0 0], [],[],[],-10);
    for boxsize_y = [14]
        max_height = boxsize_y - 1;
        tooth_heights = unique(round(linspace(max_height, 1, min(4, max_height))), 'stable');
        for tooth_height = 4       %tooth_heights %boxsize_y-1:2:1
            for removal_lenght_y = [1]
                for removal_space_y = 2*removal_lenght_y % removal_space_y= removal_lenght_y+1:2:removal_lenght_y+6
                    clearvars ct cst
                    [ct,cst] = Box_modulator(ct_no_mod,cst_no_mod, [boxsize_y 80 50], 80, [removal_lenght_y 1],[removal_space_y 30],tooth_height,132);

                    % create pln
                    for doseGridRes = [3]
                        pln = specific_pln_maker(ct, cst, doseGridRes);
                        
                        % stf
                        stf = calc_STF_Single_Energy_Layer_Beam(ct, cst, pln);
                        
                        % Calculation without optimization
                        w = ones(stf.totalNumOfBixels,1);
                        resultGUI = matRad_calcDoseForward(ct,cst,stf,pln,w);
                        resultGUI_no_mod = matRad_calcDoseForward(ct_no_mod,cst_no_mod, stf, pln, w);
                        
                        % Save with structure
                        baseFolder = 'C:\Users\victo\OneDrive\Работен плот\dskalo1\Heidelberg Physik\bachelor_thesis\dkfz\MymatRad\Project\orderly_study_26_04_data';

                        % % Build nested subfolder paths using all for-loop variables
                        % ptvFolder       = fullfile(baseFolder, sprintf('PTV_center_%d', PTV_center));
                        % boxsizeFolder   = fullfile(ptvFolder, sprintf('boxsize_y_%d', boxsize_y));
                        % toothFolder     = fullfile(boxsizeFolder, sprintf('tooth_height_%d', tooth_height));
                        % removalFolder   = fullfile(toothFolder, sprintf('removal_lenght_y_%d', removal_lenght_y));
                        % spacingFolder   = fullfile(removalFolder, sprintf('removal_space_y_%d', removal_space_y));
                        % resFolder       = fullfile(spacingFolder, sprintf('doseGridRes_%d', doseGridRes));
                        % 
                        % % Create the full path if it doesn't exist
                        % if ~exist(ptvFolder, 'dir')
                        %     mkdir(ptvFolder);
                        % end
                        % 
                        % % Create folders if they don't exist
                        % if ~exist(boxsizeFolder, 'dir')
                        %     mkdir(boxsizeFolder);
                        % end
                        % % Create folders if they don't exist
                        % if ~exist(toothFolder, 'dir')
                        %     mkdir(toothFolder);
                        % end
                        % % Create folders if they don't exist
                        % if ~exist(removalFolder, 'dir')
                        %     mkdir(removalFolder);
                        % end
                        % % Create folders if they don't exist
                        % if ~exist(spacingFolder, 'dir')
                        %     mkdir(spacingFolder);
                        % end
                        % % Create folders if they don't exist
                        % if ~exist(resFolder, 'dir')
                        %     mkdir(resFolder);
                        % end

                        % Title and analysis
                        title_text = sprintf('PTV-Center = %d; boxsize-y = %d; tooth-height = %d; removal-lenght-y = %d;\n removal-space-y = %d doseGridRes = %d HU = 132', ...
                            PTV_center, boxsize_y, tooth_height, removal_lenght_y, removal_space_y, doseGridRes);
                        data_analysis_one_Dose(ct, cst, resultGUI, pln, title_text, 50, "on");

                        % Save figure inside the deepest folder
                        % filename = sprintf('th_%d_len_%d_sp_%d_res_%d.png', tooth_height, removal_lenght_y, removal_space_y, doseGridRes);
                        % saveas(gcf, fullfile(resFolder, filename));
                        
                        % currentState = struct('PTV_center', PTV_center, ...
                        %                       'boxsize_y', boxsize_y, ...
                        %                       'tooth_height', tooth_height, ...
                        %                       'removal_lenght_y', removal_lenght_y, ...
                        %                       'removal_space_y', removal_space_y, ...
                        %                       'doseGridRes', doseGridRes);
                        % save('resume_state.mat', 'currentState');

                    end
                end
            end
        end
    end
end



%%
%matRadGUI

% To fix: boxSize_y doesn't work with odd numbers
% Make also, that the names of every file, meaning  title_text = sprintf('Cone-width=3cm Cone-spacing=5cm \n s-water-to-phantom = 35cm HU = 132');