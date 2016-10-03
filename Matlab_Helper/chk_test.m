addpath('privates');
dir_path = '~/Develop/density/test_input/output/';

% hFig = figure(2);
% clf;
% set(hFig, 'Position', [100 100 1200 900])
% Draw1stable([dir_path 'berlinvert.txt'],[dir_path 'berlinedge.txt'],'b', 1);
% 
% 
% figure;
% subplot(2,2,1);
% Draw1stable([dir_path 'testvert_full_2.txt'],[dir_path 'testedge_full_2.txt'],'r', 2);
% Draw1stable([dir_path 'testvert_thd_2.txt'],[dir_path 'testedge_thd_2.txt'],'b', 1);
% subplot(2,2,2);
% Draw1stable([dir_path 'testvert_full_5.txt'],[dir_path 'testedge_full_5.txt'],'r', 2);
% Draw1stable([dir_path 'testvert_thd_5.txt'],[dir_path 'testedge_thd_5.txt'],'b', 1);
% subplot(2,2,3);
% Draw1stable([dir_path 'testvert_full_7.txt'],[dir_path 'testedge_full_7.txt'],'r',2);
% Draw1stable([dir_path 'testvert_thd_7.txt'],[dir_path 'testedge_thd_7.txt'],'b', 1);
% subplot(2,2,4);
% Draw1stable([dir_path 'testvert_full_8.txt'],[dir_path 'testedge_full_8.txt'],'r',2);
% Draw1stable([dir_path 'testvert_thd_8.txt'],[dir_path 'testedge_thd_8.txt'],'b', 1);


% subplot(2,2,1);
% Draw1stable([dir_path 'testvert_sc.txt'],[dir_path 'testedge_sc.txt'],'b', 1);
% subplot(2,2,2);
% Draw1stable([dir_path 'testvert_mc.txt'],[dir_path 'testedge_mc.txt'],'b', 1);
% subplot(2,2,3);
% Draw1stable([dir_path 'testvert_sc_filled.txt'],[dir_path 'testedge_sc_filled.txt'],'b', 1);
% subplot(2,2,4);
% Draw1stable([dir_path 'testvert_mc_filled.txt'],[dir_path 'testedge_mc_filled.txt'],'b', 1);


%% Compare MS with disperse
dir_path = '~/Develop/density/Results/Sousbie_vs_DM/output/';
hFig = figure(1);
clf;
set(hFig, 'Position', [100 100 1200 800])
subplot(2,2,1);
DrawVTK([dir_path 'Olfactory_OP_1.fits_c10.up.NDskl.a.vtk']);
Draw1stable([dir_path 'outvert1.txt'],[dir_path 'outedge1.txt'],'b', 1);
subplot(2,2,2);
DrawVTK([dir_path 'Olfactory_OP_7.fits_c20.up.NDskl.a.vtk']);
Draw1stable([dir_path 'outvert7.txt'],[dir_path 'outedge7.txt'],'b', 1);
subplot(2,2,3);
DrawVTK([dir_path 'Olfactory_OP_8.fits_c80.up.NDskl.a.vtk']);
Draw1stable([dir_path 'outvert8.txt'],[dir_path 'outedge8.txt'],'b', 1);
subplot(2,2,4);
DrawVTK([dir_path 'Olfactory_OP_3.fits_c30.up.NDskl.a.vtk']);
Draw1stable([dir_path 'outvert3.txt'],[dir_path 'outedge3.txt'],'b', 1);