clear all;clc;close all;
a = 11.3/2;
neff = 1.447;
nclad = 1.4440;
ncore = 1.4513;DELTA=(ncore-nclad)/ncore*100;%1.4513
lambda = 1.55;
LP_type = [2,1];
m=LP_type(1);l=LP_type(2);

[E_comp, H_comp, r_s, fi_s] = ModeSolver(lambda, [a], 125/2, [1], LP_type, [0], [ncore nclad], 'LP');

PlotMode(E_comp, H_comp, r_s, fi_s, 125/2)