function [G Glevel]=gausstest2(latitude,longitude,Cn,latpdf, longpdf, T, P);

% instructions used to calculate 
% Gauss kernel. Assume no complications with the data, so that no
% significant special cases are found with lat long pairs (i.e., no points
% across the meridian of greenwich or across the equator.

% First start by assigning the grid of evaluation points.


gauss=zeros(size(T));

[mT nT]=size(T);

T=T(:);
P=P(:);
gauss=gauss(:);

% Now calculate the distances between each observation point and the data.
% Note that it is incorrect to use the latitude-longitude pairs as x-y
% couples directly, since we're really interested in distances on km. For
% this reason, the "ellipsoid vector" of an spherical Earth is introduced
% on the calculation of distances as a default. Other values depend on the
% global variable wkngplanet that is defined by runing the instruction:
% workingplanet('planetname');

global wkngplanet

if isempty(wkngplanet) | strcmp(wkngplanet,'earth')
    ellipsoid=[6371 0]
elseif strcmp(wkngplanet,'mars')
    ellipsoid=[3390 0]
elseif strcmp(wkngplanet,'mercury')
    ellipsoid=[2439 0]
elseif strcmp(wkngplanet,'moon')
    ellipsoid=[1738 0]
elseif strcmp(wkngplanet,'venus')
    ellipsoid=[6051 0]
elseif strcmp(wkngplanet,'io')
    ellipsoid=[1822 0]
end
    
    
    
for i=1:length(T),
    separation=distance([T(i) P(i)],[latitude, longitude],ellipsoid);
    gauss(i)=(1/(2*pi*Cn^2))*sum(exp((-1/(2*Cn^2))*separation.^2));
end


G=reshape(gauss, mT,nT);

Gmax=max(max(G));

Glevel=[.10*Gmax .20*Gmax .3*Gmax .4*Gmax .5*Gmax .6*Gmax .7*Gmax .8*Gmax .9*Gmax];

