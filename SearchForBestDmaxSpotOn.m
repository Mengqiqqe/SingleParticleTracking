%   SpotOn_main.m  VERSION: 1.05  Please see https://gitlab.com/tjian-darzacq-lab/spot-on-matlab for documentation
%   Copyright (C) 2017 Anders Sejr Hansen & Maxime Woringer. 
clear; clc; clearvars -global; close all; 

D_in = [0.1:0.1:1];
Dout_vs_Dmax=[];

for n = 1: numel(D_in)  
%%%%%%%%%%%%%%%%%%%%%%%% GNU LICENSE OVERVIEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% If you modify this Program, or any covered work, by linking or combining it
% with Matlab or any Matlab toolbox, the licensors of this Program grant you 
% additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, please see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% DEFINE PARAMETERS FOR MODEL-FITTING %%%%%%%%%%%%%%%%%%%%

%%%%% Choose DataSet
DataSet = 1;    % Use DataSet=1 for an example of how to process data for multiple cells from a single replicate
                % Use DataSet=2 for an example of how to process data for multiple cells from multiple replicates
data_struct = struct([]);

%%%%% Acquisition Parameters: 
type = 'ImagesD2_190ms'; % tpye of the data set (Mengqi added)
DmaxMTT = D_in(n); % the Dmax uded in MTT ttacking (Mengqi added)
TimeGap = 190; % delay between frames in milliseconds
GapsAllowed = 1; % The number of allowed gaps in the tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dZ = 0.700; % The axial observation slice in micrometers; Rougly 0.7 um for the example data (HiLo)

%%%%% Data Processing Parameters:
TimePoints = 2; % How many delays to consider: N timepoints yield N-1 delays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BinWidth = 0.010; % Bin Width for computing histogram in micrometers (only for PDF; Spot-On uses 1 nm bins for CDF)
UseEntireTraj = 1; % If UseEntireTraj=1, all dispplacements from all trajectories will be used; If UseEntireTraj=0, only the first X displacements will be used. NB. this variable was previously called UseAllTraj but has been renamed UseEntireTraj
JumpsToConsider = 4; % If UseEntireTraj=0, the first JumpsToConsiders displacements for each dT where possible will be used. 
MaxJumpPlotPDF = 3; %1.05; % the cut-off for displaying the displacement histograms plots
MaxJumpPlotCDF = 3; %3.05; % the cut-off for displaying the displacement CDF plots
MaxJump = 5; %20.05; % the overall maximal displacements to consider in micrometers
SavePlot = 1; % if SavePlot=1, key output plots will be saved to the folder "SavedPlots"; Otherwise set SavePlot = 0;
DoPlots = 1; % if DoPlots=1, Spot-On will output plots, but not if it's zero. Avoiding plots speeds up Spot-On for batch analysis

%%%%% Model Fitting Parameters:
D_Free_2State = [0.01 D_in(n)]; % min/max Diffusion constant for Free state in 2-state model (units um^2/s)
D_Bound_2State = [0.0001 0.001]; % min/max Diffusion constant for Bound state in 2-state model (units um^2/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ModelFit = 2; %Use 1 for PDF-fitting; Use 2 for CDF-fitting
DoSingleCellFit = 0; %Set to 1 if you want to analyse all single cells individually (slow). 
NumberOfStates = 2; % If NumberOfStates=2, a 2-state model will be used; If NumberOfStates=3, a 3-state model will be used 
FitIterations = 3; % Input the desired number of fitting iterations (random initial parameter guess for each)
FitLocError = 0; % If FitLocError=1, the localization error will fitted from the data
FitLocErrorRange = [0.010 0.075]; % min/max for model-fitted localization error in micrometers.
LocError = 0.035; % If FitLocError=0, LocError in units of micrometers will be used. 
UseWeights = 0; % If UseWeights=0, all TimePoints are given equal weights. If UseWeights=1, TimePoints are weighted according to how much data there is. E.g. 1dT will be weighted more than 5dT.
%D_Free_2State = [0.01 25]; % min/max Diffusion constant for Free state in 2-state model (units um^2/s)
%D_Bound_2State = [0.00001 0.0001]; % min/max Diffusion constant for Bound state in 2-state model (units um^2/s)
D_Free1_3State = [0.5 25]; % min/max Diffusion constant #1 for Free state in 3-state model (units um^2/s)
D_Free2_3State = [0.5 25]; % min/max Diffusion constant #2 for Free state in 3-state model (units um^2/s)
D_Bound_3State = [0.0001 0.05]; % min/max Diffusion constant for Bound state in 3-state model (units um^2/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% DEFINE DATA SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
if DataSet == 1 % Example DataSet 1: a single replicate of Halo-hCTCF at 134 Hz
    data_struct(1).path = sprintf('/Users/rosslab/Documents/Ross Lab/MTT&Spot-On analysis/Simulated traj analysis from Ryan/20210129 Ryan movie D2-D4/9.5ms/D2/Frame rate test/Search for best Dmax 190ms/Dmax=[0.1:0.1:1]/mat file Dmax=%d Gap=2/', D_in(n));
    
    concentration = 'ImagesD2_190ms_stack'; % file name format
    video = [1];
    
    for k = 1:numel(video)
    
        file = sprintf('%s_%03d_Tracked',concentration, video(k));
        data_struct(1).workspaces{k} = file;
    
    end
    
    data_struct(1).Include = [1:numel(video)];
    
    SampleName = sprintf('%s_Dmax=%d_Gap=%d_TP=%d', type, DmaxMTT, GapsAllowed,TimePoints);
    
    %For 20200212 data:
    %For buffer case:
    %data_struct(1).workspaces = {'buffer_133_Tracked','buffer_134_Tracked','buffer_135_Tracked','buffer_136_Tracked','buffer_137_Tracked','buffer_138_Tracked','buffer_139_Tracked','buffer_140_Tracked','buffer_141_Tracked','buffer_142_Tracked','buffer_143_Tracked','buffer_144_Tracked','buffer_145_Tracked','buffer_146_Tracked','buffer_147_Tracked','buffer_148_Tracked','buffer_149_Tracked','buffer_150_Tracked','buffer_151_Tracked','buffer_152_Tracked','buffer_153_Tracked'};
    %data_struct(1).Include = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21];
    
    %For 20191212 data:
    %For buffer case:
    %data_struct(1).workspaces = {'buffer_30ms_001_Tracked', 'buffer_30ms_002_Tracked','buffer_30ms_003_Tracked','buffer_30ms_004_Tracked','buffer_30ms_005_Tracked','buffer_30ms_006_Tracked','buffer_30ms_007_Tracked','buffer_30ms_008_Tracked','buffer_30ms_009_Tracked','buffer_30ms_010_Tracked','buffer_30ms_011_Tracked','buffer_30ms_012_Tracked','buffer_30ms_013_Tracked','buffer_30ms_014_Tracked','buffer_30ms_015_Tracked','buffer_30ms_016_Tracked','buffer_30ms_017_Tracked','buffer_30ms_018_Tracked','buffer_30ms_019_Tracked','buffer_30ms_020_Tracked','buffer_30ms_021_Tracked','buffer_30ms_022_Tracked','buffer_30ms_023_Tracked','buffer_30ms_024_Tracked','buffer_30ms_025_Tracked','buffer_30ms_026_Tracked','buffer_30ms_027_Tracked','buffer_30ms_028_Tracked','buffer_30ms_029_Tracked','buffer_30ms_030_Tracked','buffer_30ms_031_Tracked'};
    %data_struct(1).Include = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31];
    
    %For urea case:
    %data_struct(1).workspaces = {'urea_30ms_033_Tracked','urea_30ms_034_Tracked','urea_30ms_035_Tracked','urea_30ms_036_Tracked','urea_30ms_037_Tracked','urea_30ms_038_Tracked','urea_30ms_039_Tracked','urea_30ms_040_Tracked','urea_30ms_041_Tracked','urea_30ms_042_Tracked','urea_30ms_043_Tracked','urea_30ms_044_Tracked','urea_30ms_045_Tracked','urea_30ms_046_Tracked','urea_30ms_047_Tracked','urea_30ms_048_Tracked','urea_30ms_049_Tracked','urea_30ms_050_Tracked','urea_30ms_051_Tracked','urea_30ms_052_Tracked','urea_30ms_053_Tracked','urea_30ms_054_Tracked','urea_30ms_055_Tracked','urea_30ms_056_Tracked','urea_30ms_057_Tracked','urea_30ms_058_Tracked','urea_30ms_059_Tracked','urea_30ms_060_Tracked','urea_30ms_061_Tracked','urea_30ms_062_Tracked','urea_30ms_063_Tracked','urea_30ms_064_Tracked','urea_30ms_065_Tracked','urea_30ms_066_Tracked','urea_30ms_067_Tracked','urea_30ms_068_Tracked','urea_30ms_069_Tracked','urea_30ms_070_Tracked','urea_30ms_071_Tracked','urea_30ms_072_Tracked','urea_30ms_073_Tracked','urea_30ms_074_Tracked','urea_30ms_075_Tracked','urea_30ms_076_Tracked','urea_30ms_077_Tracked','urea_30ms_078_Tracked','urea_30ms_079_Tracked','urea_30ms_080_Tracked','urea_30ms_081_Tracked','urea_30ms_082_Tracked','urea_30ms_083_Tracked','urea_30ms_084_Tracked','urea_30ms_085_Tracked','urea_30ms_086_Tracked','urea_30ms_087_Tracked','urea_30ms_088_Tracked','urea_30ms_089_Tracked','urea_30ms_090_Tracked','urea_30ms_091_Tracked','urea_30ms_092_Tracked','urea_30ms_093_Tracked','urea_30ms_094_Tracked','urea_30ms_095_Tracked','urea_30ms_096_Tracked','urea_30ms_097_Tracked','urea_30ms_098_Tracked','urea_30ms_099_Tracked','urea_30ms_100_Tracked','urea_30ms_101_Tracked','urea_30ms_102_Tracked','urea_30ms_103_Tracked','urea_30ms_104_Tracked','urea_30ms_105_Tracked','urea_30ms_106_Tracked','urea_30ms_107_Tracked','urea_30ms_108_Tracked','urea_30ms_109_Tracked','urea_30ms_110_Tracked','urea_30ms_111_Tracked','urea_30ms_112_Tracked','urea_30ms_113_Tracked','urea_30ms_114_Tracked','urea_30ms_115_Tracked','urea_30ms_116_Tracked','urea_30ms_117_Tracked','urea_30ms_118_Tracked','urea_30ms_119_Tracked','urea_30ms_120_Tracked','urea_30ms_121_Tracked','urea_30ms_122_Tracked','urea_30ms_123_Tracked','urea_30ms_124_Tracked','urea_30ms_125_Tracked'};
    %data_struct(1).Include = [1:93];
  
    %SampleName = sprintf('2_%s_Dmax=%d_Gap=%d_TP=%d', type, DmaxMTT, GapsAllowed,TimePoints);
    
elseif DataSet == 2 % Example DataSet 2: two replicates of Halo-hCTCF at 134 Hz
    % 1st CTCF replicate:
    data_struct(1).path = [pwd, filesep, 'Data', filesep, 'CTCF_134Hz_rep1', filesep];
    data_struct(1).workspaces = {'U2OS_C32_Halo-CTCF_PA-JF646_1ms-633nm_134Hz_rep1_cell01', 'U2OS_C32_Halo-CTCF_PA-JF646_1ms-633nm_134Hz_rep1_cell02', 'U2OS_C32_Halo-CTCF_PA-JF646_1ms-633nm_134Hz_rep1_cell03', 'U2OS_C32_Halo-CTCF_PA-JF646_1ms-633nm_134Hz_rep1_cell04', 'U2OS_C32_Halo-CTCF_PA-JF646_1ms-633nm_134Hz_rep1_cell05'};
    data_struct(1).Include = [1,2,3,4,5];   
    % 3rd CTCF replicate:
    data_struct(2).path = [pwd, filesep, 'Data', filesep, 'CTCF_134Hz_rep3', filesep];
    data_struct(2).workspaces = {'U2OS_C32_Halo-CTCF_PA-JF646_1ms-633nm_134Hz_rep3_cell01', 'U2OS_C32_Halo-CTCF_PA-JF646_1ms-633nm_134Hz_rep3_cell02', 'U2OS_C32_Halo-CTCF_PA-JF646_1ms-633nm_134Hz_rep3_cell03', 'U2OS_C32_Halo-CTCF_PA-JF646_1ms-633nm_134Hz_rep3_cell04', 'U2OS_C32_Halo-CTCF_PA-JF646_1ms-633nm_134Hz_rep3_cell05'};
    data_struct(2).Include = [1,2,3,4,5];
    % name of merged dataset
    SampleName = 'U2OS C32 Halo-hCTCF; PA-JF646; ~134 Hz; two replicates';   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% SpotOn core mechanics %%%%%%%%%%%%%%%%%%%%%%%%%%%
Params = struct(); % Use Params to feed all the relevant data/parameters into the relevant functions
Params.TimeGap = TimeGap; Params.dZ = dZ; Params.GapsAllowed = GapsAllowed; Params.TimePoints = TimePoints; Params.BinWidth = BinWidth; Params.UseEntireTraj = UseEntireTraj; Params.DoPlots = DoPlots; Params.UseWeights = UseWeights;
Params.JumpsToConsider = JumpsToConsider; Params.MaxJumpPlotPDF = MaxJumpPlotPDF; Params.MaxJumpPlotCDF = MaxJumpPlotCDF; Params.MaxJump = MaxJump; Params.SavePlot = SavePlot; Params.ModelFit = ModelFit;
Params.DoSingleCellFit = DoSingleCellFit; Params.FitIterations = FitIterations; Params.FitLocError = FitLocError; Params.FitLocErrorRange = FitLocErrorRange; Params.LocError = LocError; Params.NumberOfStates = NumberOfStates;
Params.D_Free_2State = D_Free_2State; Params.D_Bound_2State = D_Bound_2State; Params.D_Free1_3State = D_Free1_3State; Params.D_Free2_3State = D_Free2_3State; Params.D_Bound_3State = D_Bound_3State;
Params.curr_dir = pwd; Params.SampleName = SampleName; Params.data_struct = data_struct;
% add the relevant paths
addpath(genpath([pwd, filesep, 'SpotOn_package', filesep])); 
display('Added local paths for Spot-on core mechanics');
[Output_struct] = SpotOn_core(Params);

 Dout_vs_Dmax = [Dout_vs_Dmax; D_in(n) Output_struct.merged_model_params(1)];
 
 clearvars -except D_in Dout_vs_Dmax;
 
end

Dout_vs_Dmax









