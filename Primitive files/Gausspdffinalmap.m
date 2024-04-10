function [latpdf, lonpdf, F, Flevel]=Gausspdffinalmap(latitude, longitude, Cn, plotoption)

% This version is similar to the no map version, except that it assumes
% that there is a map axes already defined, and that the results will plot
% on those axes.
% 
% Input variables now are limited to:
% latitude and longitude which are the coordinates of the data, Cn is the scale
% parameter, and plotoption, which can be 0 (no plot), 1 (plot with no data) or 2 (plot with data).
% The only plot
% produced here is a default plane, with no geographic references.
% Output still is formed by latpdf and longpdf which are the coordinates of the observation
% points used to generate the contours. 

if nargin == 3,
    plotoption = 2;
end


% To prevent the eventual case of data across the Greenwich meridian all
% the longitudes are to be expressed as longitude east of Greenwich:

if any(longitude <= 0);
    longaux=find(longitude <= 0);
    longitude(longaux)=360+longitude(longaux); 
end


% Now, it needs to find the portion of the planet where the data are found:

latrange=[min(latitude) max(latitude)];
longrange=[min(longitude) max(longitude)];

% and decide the size of the step to guarantee a given resolution within
% the area covered by the data:

resolution=200;
latstep=abs(latrange(2)-latrange(1))/(resolution+1);
longstep=abs(longrange(2)-longrange(1))/(resolution+1);

% Now find the latitude-longitude limits of the existing map projection
%figure(1);
h1=gcm;
arealatlim=h1.maplatlimit;
arealonglim=h1.maplonlimit;

% and set the grid of evaluation points

latpdf=linspace(arealatlim(1)-latstep,arealatlim(2)+latstep,resolution);
lonpdf=linspace(arealonglim(1)-longstep,arealonglim(2)+longstep,resolution);


[T,P]=meshgrid(latpdf,lonpdf);

[F, Flevel]=gausstest2(latitude,longitude,Cn,latpdf,lonpdf,T,P);


% Now  to plot the results, the commands are in the function plotpdf:

if plotoption ~= 0,
    contourfm(latpdf,lonpdf,F',Flevel);
    load('/users/ecanon/Documents/untitled folder/Matlab magnetic toolbox/Magnetic/General/pdfcolors.mat');
    set(gcf,'colormap',h)
    if plotoption == 2
        hold on
        plot(longitude,latitude,'w.','MarkerSize',10)
    end
end