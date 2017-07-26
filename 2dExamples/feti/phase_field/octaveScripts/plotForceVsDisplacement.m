clear all; close all; clc;
path = "../results_3/";
force = load([path  "force.dat"]);
displacement = load([path  "displacements.dat"]);

plot(displacement(:,2),force(:,2), '.-');
