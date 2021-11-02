%
% Make all figures for the manuscript
%
addpath('..')
clear
close all
%
% Sweeping across nutrient levels
%
figure(1)
clf
figureNutrientSweep
%
% Sweep across mixing rates and plot productivity:
%
figure(2)
clf
figureMixingRateSweep
%
% Table with biomasses
%
figure(4)
tableBiomasses
%
% Table with carbon uptakes
%
tableCarbonUptake

