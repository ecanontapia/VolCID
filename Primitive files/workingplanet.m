function workingplanet(planetname)

% This function declares the name of the planet as a global variable. It
% should be run before using Gausspdffinal or gausstest2 to get the correct
% horizontal distance between points whose position is given by lat,long
% pairs. If it is not run, those codes assume that the distribution is
% examined on Earth.

global wkngplanet

wkngplanet=planetname;