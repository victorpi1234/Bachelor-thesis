function [] = data_analysis_Dose_comparer(ct_1,cst_1,ct_2, cst_2, resultGUI_1,resultGUI_2, pln, modulator_type1, modulator_type2, rows, slice,pics_on_or_off, Doserange, Doserange_Difference)

if ~exist('rows','var') || isempty(rows)
    rows = 1; 
end
if ~exist('slice','var') || isempty(slice)
    slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct_1);
    slice = slice(3);
end
if ~exist('pics_on_or_off','var') || isempty(pics_on_or_off)
    pics_on_or_off = "on";
end
if ~exist('Doserange','var') || isempty(Doserange)
    Doserange = []; % Default: automatic color limits
end
if ~exist('Doserange_Difference','var') || isempty(Doserange_Difference)
    Doserange_Difference = []; % Default: automatic color limits
end

disp('Results function called!');
% %figure, matRad_plotSliceWrapper(gca,ct_mod_comb,cst_mod_comb,1,resultGUI.physicalDose - resultGUI_comb.physicalDose ,3,slice);
% %figure, matRad_plotSliceWrapper(gca,ct_no_mod,cst_no_mod,1,resultGUI_no_mod.physicalDose,3,slice);
% figure, matRad_plotSliceWrapper(gca,ct_mod,cst_mod,1,resultGUI.physicalDose,3,slice);
% clim([-5e-3 5e-3]); % Dose axis
% ylabel(colorbar, 'Dose (Gray)'); % Add a label to the colorbar
% title('No modulator') %With no modulator res=1
% %font size change
% fontsize = 20;
% set(gca, 'FontSize', fontsize); % Set font size for axes text (ticks, labels, etc.)
% set(get(gca, 'Title'), 'FontSize', fontsize); % Set font size for the title
% set(get(gca, 'XLabel'), 'FontSize', fontsize); % Set font size for the x-axis label
% set(get(gca, 'YLabel'), 'FontSize', fontsize); % Set font size for the y-axis label
% set(cb, 'FontSize', fontsize); % Set font size for the colorbar text
% set(get(cb, 'Label'), 'FontSize', fontsize); % Set font size for the colorbar label
%% Subplots
% Subplot 1: No comb - comb modulator
% slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct_mod);
% slice = slice(3);
figure
set(gcf, 'Visible', pics_on_or_off); 
% Subplot 2: Dose distribution with another modification (example)
subplot(rows,3,1);
matRad_plotSliceWrapper(gca, ct_1, cst_1, 1, resultGUI_1.physicalDose, 3, slice, [], [], [], ...     % thresh, alpha, contourColorMap
    [], Doserange, [], [], 'Dose [Gy]', false); 
% matRad_plotSliceWrapper(gca,ct_no_mod,cst_no_mod,1,resultGUI_no_mod.physicalDose,3,slice);
% clim([0 5e-3]);
% ylabel(colorbar, 'Dose [Gy]');
title(['Dose distribution with ' modulator_type1 ' modulator']);

subplot(rows,3,2);
matRad_plotSliceWrapper(gca, ct_2, cst_2, 1, resultGUI_2.physicalDose, 3, slice, [], [], [], ...     % thresh, alpha, contourColorMap
    [], Doserange, [], [], 'Dose [Gy]', false); 
% matRad_plotSliceWrapper(gca,ct_mod,cst_mod,1,resultGUI.physicalDose,3,slice);
%clim([0 5e-3]);
% ylabel(colorbar, 'Dose [Gy]');
title(['Dose distribution with ' modulator_type2 ' modulator']);

% Subplot 3: Another variation (example)
subplot(rows,3,3);
matRad_plotSliceWrapper(gca, ct_1, cst_1, 1, resultGUI_2.physicalDose- resultGUI_1.physicalDose, 3, slice, [], [], [], ...     % thresh, alpha, contourColorMap
    [], Doserange_Difference, [], [], 'Dose [Gy]', false); 
% matRad_plotSliceWrapper(gca,ct_1,cst_1,1,resultGUI_1.physicalDose- resultGUI_2.physicalDose ,3,slice); % Adjusted dose
% clim([-5e-3 5e-3]);
ylabel(colorbar, 'Dose [Gy]');
title('Dose difference distribution ');

% % Subplot 4: Another variation (example)
% subplot(2,2,4);
% matRad_plotSliceWrapper(gca,ct_mod,cst_mod,1,resultGUI.physicalDose - resultGUI_no_mod.physicalDose,3,slice); % Adjusted dose
% clim([-5e-3 5e-3]);
% ylabel(colorbar, 'Dose [Gy]');
% title('Comb - no modulator');
end