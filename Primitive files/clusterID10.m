function [stdofgroup, datagroupindex]=clusterID10(inputdata)

% [stdofgroup, datagroupindex]=clusterID10(inputdata)
% This function identifies clusters of geographical entities. 
% inputdata is an nx2 matrix. It MUST be in the order [LATITUDE,LONGITUDE];
% If it has longitude and latitude inverted, the results are not correct
% but the code does not generate an error.
% Each group selection has to be approved before looking for the next group, based on
% examination of the pdf contours of the prospect group.
% This provides the opportunity to use a smaller or larger standard
% deviation to get one group before moving forward. The output has two
% vectors: datagroupindex and stdofgroups.
% datagroupindex is a vector with the same number of rows as the input, 
% each containing the number of the group assigned to that element.
% stdofgroup is a vector with as many rows as groups identified, each row
% having the std that was used to accept the group graphically.


% Test for duplicate entries on the data, and if found, then delete them

datatemp=unique(inputdata,'rows');

if length(datatemp) ~= length(inputdata)
    warning('Duplicate locations were found and removed')
    inputdata=datatemp;
end

% To prevent the eventual case of data across the Greenwich meridian all
% the longitudes are to be expressed as longitude east of Greenwich:

if any(inputdata(:,2) <= 0);
    longaux=find(inputdata(:,2) <= 0);
    inputdata(longaux,2)=360+inputdata(longaux,2); 
end

% As a default, the map figure will include a margin around the data that
% is equal to 10% of the width of the space used by the data.

latlim=[min(inputdata(:,1)), max(inputdata(:,1))];
lonlim=[min(inputdata(:,2)), max(inputdata(:,2))];
map_lat_extension=latlim(2)-latlim(1);
map_lon_extension=lonlim(2)-lonlim(1);
map_lat_lim=[latlim(1)-0.1*map_lat_extension, latlim(2)+0.1*map_lat_extension];
map_lon_lim=[lonlim(1)-0.1*map_lon_extension, lonlim(2)+0.1*map_lon_extension];



% Separation is calculated for Earth as default, but the radius can be adjusted to 
% reflect other planets by running workingplanet on the command window first,
% where the planet is explicitly stated. 

global wkngplanet

if isempty(wkngplanet) | strcmp(wkngplanet,'earth')
    ellipsoid=[6371 0];
elseif strcmp(wkngplanet,'mars')
    ellipsoid=[3390 0];
elseif strcmp(wkngplanet,'mercury')
    ellipsoid=[2439 0];
elseif strcmp(wkngplanet,'moon')
    ellipsoid=[1738 0];
elseif strcmp(wkngplanet,'venus')
    ellipsoid=[6051 0];
elseif strcmp(wkngplanet,'io')
    ellipsoid=[1822 0];
end

% initialize some variables

datagroupindex=zeros(length(inputdata),1);
groupnumber=0;
flowcontrol1=0;
stdofref=1;

% First, calculate the distances between points
% The ith column of separation has the distances of the ith observation to 
% all the observations on the database, including itself.
[m,n]=size(inputdata);

separation=zeros(m,m);

for i=1:m
    separation(:,i)=distance(inputdata(i,1),inputdata(i,2),inputdata(:,1),inputdata(:,2),ellipsoid);
end

separation(separation == 0)=NaN;

% Now start the search for groups

while flowcontrol1 == 0

    groupnumber=groupnumber+1;
    
    % First select the location of the seed. 

  [x,y]=min(mean(separation,'omitnan'));
    

% Remember: y is the index of the seed location for the group in turn


% Now the nearest neighbour distribution is obtained

    NND=min(separation);
    
    
     f1=figure(1);
            if groupnumber == 1
                set(f1,'Position', [999   688   978   651])  %Mac desk
            end
            counts=hist(log(NND));
            hhist0=histogram(log(NND));
            hold on
            xstd=[mean(log(NND),'omitnan')
            mean(log(NND),'omitnan') + std(log(NND),'omitnan')
            mean(log(NND),'omitnan') + 1.5*std(log(NND),'omitnan')
            mean(log(NND),'omitnan') + 1.75*std(log(NND),'omitnan')
            mean(log(NND),'omitnan') + 2*std(log(NND),'omitnan')];
            xstd=[xstd xstd];
            ystd=[0 0.99*max(counts)
            0 0.95*max(counts)
            0 0.90*max(counts)
            0 0.85*max(counts)
            0 0.80*max(counts)];
             
            for i=1:length(ystd)
                plot(xstd(i,:),ystd(i,:),'r-')
            end
            text(xstd(1), ystd(1,2),'Mean')
            text(xstd(2), ystd(2,2),'1 std')
            text(xstd(3), ystd(3,2),'1.5 std')
            text(xstd(4), ystd(4,2),'1.75 std')
            text(xstd(5), ystd(5,2),'2 std')
            
            xlabel('Log NND')


% The default separation distance of one cluster is taken as the
% average NND plus 1 standard deviations. This means that more than 68% of the NND
% are taken into consideration. Note that two std = 95% and 3 std = 99.73%

  f2=figure(2);

  
  if groupnumber == 1
      set(f2,'Position',[1223  12   978    651]) %Mac desk
      hmap=axesm('mercator','MapLatLim',map_lat_lim,'MapLonLim',map_lon_lim,'MeridianLabel','on','ParallelLabel','on','Frame','on');
      setm(hmap,'MLabelLocation',map_lon_lim,'MLabelParallel','south','PLabelLocation',map_lat_lim,'PLabelMeridian','west','MLabelRound',-2,'PLabelRound',-2)
      plotm(inputdata(:,1),inputdata(:,2),'r^','MarkerFaceColor','r','MarkerSize',4)
      
      hold on
  end
  plotm(inputdata(y,1),inputdata(y,2),'bo')
   

flowcontrolaux=0; % Here starts the main loop

while flowcontrolaux == 0;
    clusterdistance=exp(mean(log(NND),'omitnan')+ stdofref*std(log(NND),'omitnan'));

% Start defining one cluster by identifying all the data that are at a
% distance equal or less than the selected typical cluster distance from
% the "center", and gradually increase the search to all the new added
% data. The search stops when the distance from one point in the cluster
% to another not in the cluster exceeds the typical scale defined above.


    datagroupindex(y)=groupnumber;
   

    flowcontrol=0;
    iterationnumber=1;
    newdata=y;
    separation1=separation;

    while flowcontrol == 0 % This loop extends one group around a seed
  
 
        for i=1:length(newdata)
            t1=find(separation1(newdata(i),:) <= clusterdistance);
            if i == 1
                newdataaux=t1;
            else
                newdataaux=[newdataaux, t1];
            end
            separation1(newdata(i),:)=NaN;
            separation1(:,newdata(i))=NaN;
        end
        newdata=unique(newdataaux);
        datagroupindex(newdata)=groupnumber;

        if isempty(newdata)
            flowcontrol = 1;
        end

    end  % This is the end of the while loop that extends one group around a given seed
    
 
            grouptobeplotted=inputdata(datagroupindex==groupnumber,:);
            [m1,n1]=size(grouptobeplotted);
            separationaux=zeros(m1,m1);
            for j1=1:m1
                separationaux(:,j1)=distance(grouptobeplotted(j1,1),grouptobeplotted(j1,2),grouptobeplotted(:,1),grouptobeplotted(:,2),ellipsoid);
            end
            separationaux(separationaux == 0)=NaN;
            NNDaux=min(separationaux);
            f1=figure(1);
            if groupnumber == 1
                set(f1,'Position', [999   688   978   651])  %Mac desk
            end
            if length(grouptobeplotted) <= 2
                warning('Too few points in a group. No PDF will be calculated')
                handlehistaux=[];
                counts=[];
                handlehistlines=[];
                handlehisttext=[];
            else
                counts=hist(log(NNDaux));
                handlehistaux=histogram(log(NNDaux));
                set(handlehistaux,'FaceAlpha',.5,'FaceColor','r')
                handlehistaux.BinWidth=hhist0.BinWidth;
                handlehistaux.NumBins=hhist0.NumBins;
                handlehistaux.BinEdges=hhist0.BinEdges;
                xstd=[mean(log(NNDaux),'omitnan')
                    mean(log(NNDaux),'omitnan')
                    NaN
                mean(log(NNDaux),'omitnan') + std(log(NNDaux),'omitnan')
                mean(log(NNDaux),'omitnan') + std(log(NNDaux),'omitnan')
                NaN
                mean(log(NNDaux),'omitnan') + 1.5*std(log(NNDaux),'omitnan')
                mean(log(NNDaux),'omitnan') + 1.5*std(log(NNDaux),'omitnan')
                NaN
                mean(log(NNDaux),'omitnan') + 1.75*std(log(NNDaux),'omitnan')
                mean(log(NNDaux),'omitnan') + 1.75*std(log(NNDaux),'omitnan')
                NaN
                mean(log(NNDaux),'omitnan') + 2*std(log(NNDaux),'omitnan')
                mean(log(NNDaux),'omitnan') + 2*std(log(NNDaux),'omitnan')];
                
                ystd=[0
                    0.99*max(counts)
                    NaN 
                0
                0.95*max(counts)
                NaN 
                0
                0.90*max(counts)
                NaN
                0
                0.85*max(counts)
                NaN 
                0
                0.80*max(counts)];
                handlehistlines=plot(xstd,ystd,'y-','LineWidth',4);
                
                xstdtext=[xstd(1), xstd(4), xstd(7), xstd(10), xstd(13)];
                ystdtext=[ystd(2), ystd(5), ystd(8), ystd(11), ystd(14)];
                tstdtext={'Mean', '1 std', '1.5 std', '1.75 std', '2 std'};
                handlehisttext=text(xstdtext, ystdtext,tstdtext,'Color','y');
                
                
            
                info1=['N of the group = ',num2str(length(grouptobeplotted))];
                title(info1);
                xlabel('Log NND')
                hold off
            end
            
            figure(2)
            grouphandle=plotm(inputdata(datagroupindex==groupnumber,1),inputdata(datagroupindex==groupnumber,2),'kv','MarkerFaceColor',rand(1,3),'MarkerSize',15);
            drawnow
            if length(grouptobeplotted) < 10
                warning('Too few points in a group. No PDF will be calculated')
                htest=[];
                c=[];
            else
                Cn=max(mean(separationaux,'omitnan'))/6;
                [latpdf lonpdf F Flevel]=Gausspdffinalmap(grouptobeplotted(:,1),grouptobeplotted(:,2),Cn,0);
                [c,htest]=contourm(latpdf,lonpdf,F',Flevel,'Linewidth',1);
                drawnow
                
            end
            
      
    group_ok=input('Is the group OK? [y], [n]   ','s');
        if strcmp(group_ok,'n')
            delete(grouphandle)
            delete(htest)
            delete(handlehistaux)
            delete(handlehistlines)
            delete(handlehisttext)
            clear counts
            clear grouptobeplotted
            clear separationaux
            clear c
            [ma na]=find(datagroupindex == groupnumber);
            datagroupindex(ma)=groupnumber-1;
            stdofref=input('Enter the new value of std ');
            figure(1)
            hold on
        else
            delete(htest)
            clear grouptobeplotted
            clear separationaux
            clear c
            separation=separation1;
            flowcontrolaux=1;
            if groupnumber == 1
                stdofgroup=stdofref;
            else
                stdofgroup=[stdofgroup; stdofref];
             
            end
            stdofref=1;
        end
end
        
    if all(all(isnan(separation)))
        flowcontrol1 = 1;
    end

end




    

