clear all; close all; clc;

disp = load("../results_strong_scaling1/displacements.dat");
force = load("../results_strong_scaling1/force.dat");

plot(-disp(:,2), -force(:,2));