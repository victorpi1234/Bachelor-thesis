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
clearvars ct cst
for PTV_center = [0]
    [ct_no_mod, cst_no_mod] = makeCST_water(ct_no_mod, 30, 20, [PTV_center 0 0], [],[],[],-80);
    for tooth_height = 11 %1:15
        for boxsize_y =  tooth_height + 5
            for removal_lenght_y = 1 %[1, 2]
                for removal_space_y = 2*removal_lenght_y % removal_space_y= removal_lenght_y+1:2:removal_lenght_y+6
                    clearvars ct cst
                    [ct,cst] = Box_modulator(ct_no_mod,cst_no_mod, [boxsize_y 80 50], 90, [removal_lenght_y 1],[removal_space_y 30],tooth_height,132);

                    % create pln
                    for doseGridRes = [1]
                        pln = specific_pln_maker(ct, cst, doseGridRes);
                        
                        % stf
                        stf = calc_STF_Single_Energy_Layer_Beam(ct, cst, pln);
                        
                        % Calculation without optimization
                        w = ones(stf.totalNumOfBixels,1);
                        resultGUI = matRad_calcDoseForward(ct,cst,stf,pln,w);
                        resultGUI_no_mod = matRad_calcDoseForward(ct_no_mod,cst_no_mod, stf, pln, w);
                        
                        % % Save with structure
                        % baseFolder = 'C:\Users\victo\OneDrive\Работен плот\dskalo1\Heidelberg Physik\bachelor_thesis\dkfz\MymatRad\Project\orderly_study_26_04_data';
                        % 
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

                        % Let's plot single profiles that are perpendicular to the beam direction
                        ixProfileX = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
                        ixProfileX = ixProfileX(1)-8:ixProfileX(1)+8;
                        ixProfileY = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
                        for ixProfileY = ixProfileY(2):ixProfileY(2)+1
                            slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
                            slice = slice(3);
                            profile = resultGUI.physicalDose(:,ixProfileY,slice);
                            % Store both profiles for comparison
                            if ixProfileY == ixProfileY(2)
                                profile1 = profile;
                                [max1, idx1] = max(profile1);
                            else
                                profile2 = profile;
                                [max2, idx2] = max(profile2);
        
                            % Calculate and display the shift
                                shift = abs(idx2 - idx1);
                                disp(['Peak shift between profiles: ' num2str(shift) ' voxels']);
                            end

                            figure,plot(profile,'LineWidth',2),grid on,hold on,
                            plot(profile,'LineWidth',2),legend({'original profile'}),
                            xlabel('mm'),ylabel('Gy(RBE)'),title('profile plot')
                        end

                        % % Title and analysis
                        % title_text = sprintf('PTV-Center = %d; boxsize-y = %d; tooth-height = %d; removal-lenght-y = %d;\n removal-space-y = %d doseGridRes = %d HU = 132', ...
                        %     PTV_center, boxsize_y, tooth_height, removal_lenght_y, removal_space_y, doseGridRes);
                        % data_analysis_one_Dose(ct, cst, resultGUI, pln, title_text, 50);
                        % 
                        % % Save figure inside the deepest folder
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
% Let's compare the new recalculation against the optimization result.
plane = 3;
doseWindow = [0 max([resultGUI.RBExDose(:); resultGUI_isoShift.RBExDose(:)])];

figure,title('original plan')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('shifted plan')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_isoShift.RBExDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);

absDiffCube = resultGUI.RBExDose-resultGUI_isoShift.RBExDose;
figure,title('absolute difference')
matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);

% Let's plot single profiles that are perpendicular to the beam direction
ixProfileY = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
ixProfileY = ixProfileY(2);
profileOrginal = resultGUI.RBExDose(:,ixProfileY,slice);
profileShifted = resultGUI_isoShift.RBExDose(:,ixProfileY,slice);

figure,plot(profileOrginal,'LineWidth',2),grid on,hold on,
plot(profileShifted,'LineWidth',2),legend({'original profile','shifted profile'}),
xlabel('mm'),ylabel('Gy(RBE)'),title('profile plot')
%%
matRadGUI
