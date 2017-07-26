close all; clear all; clc;

force = load("../results_wedge_split_test/forceTop.dat");
displacements = load("../results_wedge_split_test/displacementsTop.dat");

plot(displacements(:,2), force(:,2),'.-')
hold on;

force = load("../results_wedge_split_test/forceBottom.dat");
displacements = load("../results_wedge_split_test/displacementsBottom.dat");

plot(displacements(:,2), force(:,2),'.-')

xlabel("displacements");
ylabel("force");