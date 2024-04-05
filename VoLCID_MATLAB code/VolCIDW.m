function varargout = VolCIDW(varargin)
% VOLCIDW MATLAB code for VolCIDW.fig
%      VOLCIDW, by itself, creates a new VOLCIDW or raises the existing
%      singleton*.
%
%      H = VOLCIDW returns the handle to a new VOLCIDW or the handle to
%      the existing singleton*.
%
%      VOLCIDW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOLCIDW.M with the given input arguments.
%
%      VOLCIDW('Property','Value',...) creates a new VOLCIDW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VolCIDW_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VolCIDW_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VolCIDW

% Last Modified by GUIDE v2.5 09-Jul-2020 11:38:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VolCIDW_OpeningFcn, ...
                   'gui_OutputFcn',  @VolCIDW_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before VolCIDW is made visible.
function VolCIDW_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VolCIDW (see VARARGIN)

set([handles.pushbutton_ChMembership, handles.pushbutton_NextGroup, handles.pushbuttonStartGrouping, handles.pushbuttonFinishGroupSearch],'Enable','off')

% Choose default command line output for VolCIDW
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes VolCIDW wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VolCIDW_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1



% --- Executes on button press in SelectFile.
function SelectFile_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global inputdata_cell 

[handles.filename, pathname] = uigetfile('*.*');

shflnm=length(handles.filename);

if strcmp(handles.filename(shflnm),'x') || strcmp(handles.filename(shflnm),'s')
    handles.inputdata=xlsread(handles.filename);
elseif strcmp(handles.filename(shflnm-2:shflnm),'mat')
    inputdata_cell=load(fullfile(pathname, handles.filename));
    inputdata_fieldnames=fieldnames(inputdata_cell);
    t1=length(inputdata_fieldnames);
    in_names=inputdata_fieldnames;
    in_names{t1+1,1}='Select Data';
    in_names=flipud(in_names);
    set(handles.VariableNamesMenu, 'string', in_names)
    set(handles.popupmenu_BaseImage, 'string', in_names)
    set(handles.popupmenu_Gref, 'string', in_names)

    if length(inputdata_fieldnames) == 1
        handles.text3.String=[];
        inputdata_varname=inputdata_fieldnames{1};
        handles.inputdata=inputdata_cell.(inputdata_varname);
    else
        handles.inputdata=[];
    end
    
elseif strcmp(handles.filename(shflnm-2:shflnm),'csv')
    handles.inputdata=csvread(handles.filename);
elseif strcmp(handles.filename(shflnm-2:shflnm),'txt')
    handles.inputdata=textread(handles.filename, '%f %f');
end

if isempty(handles.inputdata) == 0
    handles.inputdata=CreateAxes(handles.inputdata, handles.axes1);
    handles.alldata=plotm(handles.inputdata(:,1),handles.inputdata(:,2),'r^','MarkerFaceColor','r','MarkerSize',4);
    handles.hhist0=mkhistogram(handles.inputdata,handles.axes2);
    total_obs=num2str(length(handles.inputdata));
    total_obs_info=['Number of Observations: ',total_obs];
    set(handles.text_TotalObsInfo,'string',total_obs_info)
    handles.histlinear0=mklinhistogram(handles.inputdata, handles.axes3);
end
    
    
 guidata(hObject,handles)



% --- Executes on selection change in VariableNamesMenu.
function VariableNamesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to VariableNamesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns VariableNamesMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from VariableNamesMenu

global inputdata_cell 

val=get(hObject,'Value');
inputdata_varname=handles.VariableNamesMenu.String{handles.VariableNamesMenu.Value};
handles.inputdata=inputdata_cell.(inputdata_varname);
[m,n]=size(handles.inputdata);
if (n == 2) | (m == 2)
    handles.inputdata=CreateAxes(handles.inputdata, handles.axes1);
    handles.alldata=plotm(handles.inputdata(:,1),handles.inputdata(:,2),'r^','MarkerFaceColor','r','MarkerSize',4);
    handles.text3.String='';
    handles.hhist0=mkhistogram(handles.inputdata,handles.axes2);
    total_obs=num2str(length(handles.inputdata));
    total_obs_info=['Number of Observations: ',total_obs];
    set(handles.text_TotalObsInfo,'string',total_obs_info)
    handles.histlinear0=mklinhistogram(handles.inputdata, handles.axes3);
    handles.Gref=[];
    set(handles.pushbuttonStartGrouping,'Enable','on')
else
    handles.text3.String='Need to select a different variable from the above menu';
end
guidata(hObject,handles)
    

% --- Executes during object creation, after setting all properties.
function VariableNamesMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VariableNamesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [inputdata, hhist0, histlinear0]=DefineVariableName(handles)

global inputdata_cell 

inputdata_varname=SelectVariableName(handles);
inputdata=inputdata_cell.(inputdata_varname);
[m,n]=size(inputdata);
if (n == 2) || (m == 2)
    handles.inputdata=CreateAxes(inputdata, handles.axes1);
    handles.alldata=plotm(handles.inputdata(:,1),handles.inputdata(:,2),'r^','MarkerFaceColor','r','MarkerSize',4);
    hhist0=mkhistogram(inputdata,handles.axes2);
    total_obs=num2str(length(inputdata));
    total_obs_info=['Number of Observations: ',total_obs];
    set(handles.text_TotalObsInfo,'string',total_obs_info)
    handles.histlinear0=mklinhistogram(handles.inputdata, handles.axes3);
    handles.Gref=[];
else
    hhist0=[];
    histlinear0=[];
    handles.text3.String='Need to select a different variable from the above menu';
end



function inputdata_varname=SelectVariableName(handles)

inputdata_varname=handles.VariableNamesMenu.String{handles.VariableNamesMenu.Value};


% --- Executes on selection change in popupmenuPlanet.
function popupmenuPlanet_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuPlanet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuPlanet contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuPlanet

global ellipsoid

wkngplanet=handles.popupmenuPlanet.String{handles.popupmenuPlanet.Value};
if isempty(wkngplanet) | strcmp(wkngplanet,'Earth')
    ellipsoid=[6371 0];
elseif strcmp(wkngplanet,'Mars')
    ellipsoid=[3390 0];
elseif strcmp(wkngplanet,'Mercury')
    ellipsoid=[2439 0];
elseif strcmp(wkngplanet,'Moon')
    ellipsoid=[1738 0];
elseif strcmp(wkngplanet,'Venus')
    ellipsoid=[6051 0];
elseif strcmp(wkngplanet,'Io')
    ellipsoid=[1822 0];
end


ellip_info=['Planet Radius: ',num2str(ellipsoid(1,1)), ' km'];
set(handles.text_Ellipsoid,'String',ellip_info)


% --- Executes on button press in pushbuttonStartGrouping.
function pushbuttonStartGrouping_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStartGrouping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ellipsoid  


wkngplanet=handles.popupmenuPlanet.String{handles.popupmenuPlanet.Value};
if isempty(wkngplanet) || strcmp(wkngplanet,'Earth')
    ellipsoid=[6371 0];
end

ellip_info=['Planet Radius: ',num2str(ellipsoid(1,1)), ' km'];
set(handles.text_Ellipsoid,'String',ellip_info)

handles.separation=distance_matrix(handles.inputdata);
handles.datagroupindex=zeros(length(handles.inputdata),1);  
handles.stdofref=1;
handles.groupnumber=1;

handles.y=seed_location(handles.separation, handles);

handles.NND=min(handles.separation);
[grouptobeplotted, handles.datagroupindex, handles.separationaux, handles.separation1, NNDaux, handles.grouphandle, handles.clusterdistance] =findgroup(handles.stdofref, handles.NND, handles.groupnumber, handles.y, handles);

if length(grouptobeplotted) >= 5
    info2=['Calculating group PDF, please wait'];
    set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
    handles.htest=Gausskerneltemp(grouptobeplotted, handles.separationaux, handles);
    [handles.handlehistaux, handles.histaux_lines, handles.histaux_txt]=grouphist(handles.hhist0, grouptobeplotted, NNDaux, handles);
    [handles.handlehistlinaux, handles.histlinaux_lines, handles.histlinaux_txt]=grouphistlin(handles.histlinear0, grouptobeplotted, NNDaux, handles);
    info2=['Please confirm membership or repeat with new std'];
    set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
else
    handles.htest=[];
    handles.handlehistaux=[];
    handles.histaux_lines=[];
    handles.histaux_txt=[];
    handles.handlehistlinaux=[];
    handles.histlinaux_lines=[];
    handles.histlinaux_txt=[];
    info2=['Too few points in current group. No PDF will be calculated'];
    set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
    info5=['Too few points in a group. No group histrogram shown'];
    set(handles.text_NoHistInfo,'string',info5)
end

info4=['Groups identified so far: ', num2str(handles.groupnumber)];
set(handles.text_NumberofGroupInfo,'string',info4)
info1=['N of the group = ',num2str(length(grouptobeplotted))];
set(handles.text_CurrentGroupInfo,'string',info1)
info1a=['Equivalent linear distance = ',num2str(handles.clusterdistance), '  km'];
set(handles.text_LinearDistance,'string',info1a)
non_clustered= length(find(handles.datagroupindex == 0));
info4a=['Remaining non-Grouped observations: ', num2str(non_clustered)];
set(handles.text_NumberOfNonClusterData,'string',info4a)

set(handles.pushbuttonStartGrouping,'Enable','off')
set([handles.pushbutton_ChMembership, handles.pushbutton_NextGroup, handles.pushbuttonFinishGroupSearch],'Enable','on')
guidata(hObject,handles)




function separation=distance_matrix(inputdata)

global ellipsoid

[m,n]=size(inputdata);

separation=zeros(m,m);

for i=1:m
    separation(:,i)=distance(inputdata(i,1),inputdata(i,2),inputdata(:,1),inputdata(:,2),ellipsoid);
end

separation(separation == 0)=NaN;

    

function y=seed_location(separation1, handles)

  
    
[x,y]=min(mean(separation1,'omitnan'));
axes(handles.axes1);
plotm(handles.inputdata(y,1),handles.inputdata(y,2),'bo','MarkerSize',10)
drawnow




function [grouptobeplotted, datagroupindex, separationaux, separation1, NNDaux, grouphandle, clusterdistance] =findgroup(stdofref, NND, groupnumber, y, handles)

global  ellipsoid 

clusterdistance=exp(mean(log(NND),'omitnan')+ stdofref*std(log(NND),'omitnan'));
  
handles.datagroupindex(y)=groupnumber;
flowcontrol=0;

newdata=y;
separation1=handles.separation;

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
        handles.datagroupindex(newdata)=groupnumber;

        if isempty(newdata)
            flowcontrol = 1;
        end

end 

datagroupindex=handles.datagroupindex;
grouptobeplotted=handles.inputdata(handles.datagroupindex==groupnumber,:);
[m1,n1]=size(grouptobeplotted);
separationaux=zeros(m1,m1);
        for j1=1:m1
              separationaux(:,j1)=distance(grouptobeplotted(j1,1),grouptobeplotted(j1,2),grouptobeplotted(:,1),grouptobeplotted(:,2),ellipsoid);
         end
separationaux(separationaux == 0)=NaN;
NNDaux=min(separationaux);
axes(handles.axes1)  
grouphandle=plotm(handles.inputdata(handles.datagroupindex==groupnumber,1),handles.inputdata(handles.datagroupindex==groupnumber,2),'kv','MarkerFaceColor',rand(1,3),'MarkerSize',15);
drawnow



function htest=Gausskerneltemp(grouptobeplotted, separationaux, handles)

axes(handles.axes1)

pdfcontourcolors=[.8, .8, .8;
                    0, .2, 1;
                    0, .7, 1;
                    0, 1, .7;
                    .2, .9, .1;
                    1, 1, .2;
                    1, .8, .2;
                    1, .5, .5;
                    1, .1, .2;
                    1, 0, 1];
  
            if length(grouptobeplotted) < 6
                info2=['Too few points in current group. No PDF will be calculated'];
                set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
                htest=[];
                c=[];
            else
                set(handles.text_PDFInfo,'string','')
                Cn=max(mean(separationaux,'omitnan'))/6;
                [latpdf lonpdf F Flevel]=Gausspdffinalmap(grouptobeplotted(:,1),grouptobeplotted(:,2),Cn,0);
                [c,htest]=contourm(latpdf,lonpdf,F',Flevel(2:10),'Linewidth',1);
                hpatch = get(htest,'Children');
                for j = 1:9,
                    ch = hpatch(j);
                    set(ch,'Color',pdfcontourcolors(11-j,:));
                end
    
            end
            
            

function [handlehistaux, histaux_lines, histaux_txt]=grouphist(hhist0, grouptobeplotted, NNDaux, handles) 

axes(handles.axes2)
            
         if length(grouptobeplotted) <= 5
              info5=['Too few points in a group. No histrogram shown'];
              set(handles.text_NoHistInfo,'string',info5)
              handlehistaux=[];
              histaux_lines=[];
              histaux_txt=[];
              counts=[];
              lines=[];
              
         else
             set(handles.text_NoHistInfo,'string','')
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
               
                histaux_lines=plot(xstd,ystd,'y-','LineWidth',4);
                
                xstdtext=[xstd(1), xstd(4), xstd(7), xstd(10), xstd(13)];
                ystdtext=[ystd(2), ystd(5), ystd(8), ystd(11), ystd(14)];
                tstdtext={'Mean', '1 std', '1.5 std', '1.75 std', '2 std'};
                histaux_txt=text(xstdtext, ystdtext,tstdtext,'Color','y');
                
         end
         
        
         
            
function [handlehistlinaux, histlinaux_lines, histlinaux_txt]=grouphistlin(histlinear0, grouptobeplotted, NNDaux, handles) 

axes(handles.axes3)
            
         if length(grouptobeplotted) <= 5
              handlehistlinaux=[];
              histlinaux_lines=[];
              histlinaux_txt=[];
              counts=[];
              lines=[];
              
         else
              counts=hist(NNDaux);
              handlehistlinaux=histogram(NNDaux);
              set(handlehistlinaux,'FaceAlpha',.5,'FaceColor','r')
              handlehistlinaux.BinWidth=histlinear0.BinWidth;
              handlehistlinaux.NumBins=histlinear0.NumBins;
              handlehistlinaux.BinEdges=histlinear0.BinEdges;
              xstd=[mean(NNDaux,'omitnan')
                    mean(NNDaux,'omitnan')
                    NaN
                    mean(NNDaux,'omitnan') + std(NNDaux,'omitnan')
                    mean(NNDaux,'omitnan') + std(NNDaux,'omitnan')
                    NaN
                    mean(NNDaux,'omitnan') + 1.5*std(NNDaux,'omitnan')
                    mean(NNDaux,'omitnan') + 1.5*std(NNDaux,'omitnan')
                    NaN
                    mean(NNDaux,'omitnan') + 1.75*std(NNDaux,'omitnan')
                    mean(NNDaux,'omitnan') + 1.75*std(NNDaux,'omitnan')
                    NaN
                    mean(NNDaux,'omitnan') + 2*std(NNDaux,'omitnan')
                    mean(NNDaux,'omitnan') + 2*std(NNDaux,'omitnan')];
                
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
               
                histlinaux_lines=plot(xstd,ystd,'y-','LineWidth',4);
                
                xstdtext=[xstd(1), xstd(4), xstd(7), xstd(10), xstd(13)];
                ystdtext=[ystd(2), ystd(5), ystd(8), ystd(11), ystd(14)];
                tstdtext={'Mean', '1 std', '1.5 std', '1.75 std', '2 std'};
                histlinaux_txt=text(xstdtext, ystdtext,tstdtext,'Color','y');
                
         end
         
         
            


% --- Executes on button press in pushbutton_NextGroup.
function pushbutton_NextGroup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_NextGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isempty(handles.handlehistaux) == 0
    delete(handles.histaux_txt)
    delete(handles.histaux_lines)
    delete(handles.handlehistaux)
    delete(handles.histlinaux_txt)
    delete(handles.histlinaux_lines)
    delete(handles.handlehistlinaux)
    delete(handles.htest)
end

[ma, na]=find(handles.datagroupindex == 0);

if isempty(ma)
    info2=['All observations have been grouped'];
    set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
else
    handles.separation=handles.separation1;
    if handles.groupnumber == 1
        handles.stdofgroup=handles.stdofref;
        handles.LinearDistance=handles.clusterdistance;
    else
        handles.stdofgroup=[handles.stdofgroup; handles.stdofref];
        handles.LinearDistance=[handles.LinearDistance; handles.clusterdistance];
    end

    handles.stdofref=1;
    set(handles.edit_StdofRef, 'String', num2str(handles.stdofref));
    handles.y=seed_location(handles.separation, handles);
    handles.groupnumber=handles.groupnumber+1;

    handles.NND=min(handles.separation);
    [grouptobeplotted, handles.datagroupindex, handles.separationaux, handles.separation1, NNDaux, handles.grouphandle, handles.clusterdistance] =findgroup(handles.stdofref, handles.NND, handles.groupnumber, handles.y, handles);
    if length(grouptobeplotted) >= 5
        info2=['Calculating group PDF, please wait'];
        set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
        handles.htest=Gausskerneltemp(grouptobeplotted, handles.separationaux, handles);
        [handles.handlehistaux, handles.histaux_lines, handles.histaux_txt]=grouphist(handles.hhist0, grouptobeplotted, NNDaux, handles);
        [handles.handlehistlinaux, handles.histlinaux_lines, handles.histlinaux_txt]=grouphistlin(handles.histlinear0, grouptobeplotted, NNDaux, handles);
        info2=['Please confirm membership or repeat with new std'];
        set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
    else
    
        info2=['Too few points in current group. No PDF will be calculated'];
        set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
        info5=['Too few points in a group. No histrogram shown'];
        set(handles.text_NoHistInfo,'string',info5)
    end

    info4=['Groups identified so far: ', num2str(handles.groupnumber)];
    set(handles.text_NumberofGroupInfo,'string',info4)
    info1=['N of the group = ',num2str(length(grouptobeplotted))];
    set(handles.text_CurrentGroupInfo,'string',info1)
    info1a=['Equivalent linear distance = ',num2str(handles.clusterdistance), '  km'];
    set(handles.text_LinearDistance,'string',info1a)
    non_clustered= length(find(handles.datagroupindex == 0));
    info4a=['Remaining non-Grouped observations: ', num2str(non_clustered)];
    set(handles.text_NumberOfNonClusterData,'string',info4a)
end

guidata(hObject,handles)


% --- Executes on button press in pushbutton_ChMembership.
function pushbutton_ChMembership_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ChMembership (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



if isempty(handles.handlehistaux) == 0
    delete(handles.histaux_txt)
    delete(handles.histaux_lines)
    delete(handles.handlehistaux)
    delete(handles.histlinaux_txt)
    delete(handles.histlinaux_lines)
    delete(handles.handlehistlinaux)
    delete(handles.htest)
    delete(handles.grouphandle) 
end

[ma na]=find(handles.datagroupindex == handles.groupnumber);
handles.datagroupindex(ma,1)=0;

handles.stdofref=get(handles.edit_StdofRef, 'String');
handles.stdofref=str2num(handles.stdofref);

 axes(handles.axes1);
[grouptobeplotted, handles.datagroupindex, handles.separationaux, handles.separation1, NNDaux, handles.grouphandle, handles.clusterdistance]=findgroup(handles.stdofref, handles.NND, handles.groupnumber, handles.y, handles);

if length(grouptobeplotted) >= 6
    info2=['Calculating group PDF, please wait'];
    set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
    handles.htest=Gausskerneltemp(grouptobeplotted, handles.separationaux, handles);
    [handles.handlehistaux, handles.histaux_lines, handles.histaux_txt]=grouphist(handles.hhist0, grouptobeplotted, NNDaux, handles);
    [handles.handlehistlinaux, handles.histlinaux_lines, handles.histlinaux_txt]=grouphistlin(handles.histlinear0, grouptobeplotted, NNDaux, handles);
    info2=['Please confirm membership or repeat with new std'];
    set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
else
    handles.htest=[];
    handles.handlehistaux=[];
    handles.histaux_txt=[];
    handles.histaux_lines=[];
    handles.handlehistlinaux=[];
    handles.histlinaux_lines=[];
    handles.histlinaux_txt=[];
    info2=['Too few points in current group. No PDF will be calculated'];
    set(handles.text_PDFInfo,'string',info2,'FontSize',12,'FontWeight','Bold')
    info5=['Too few points in a group. No group histrogram shown'];
    set(handles.text_NoHistInfo,'string',info5)
end

info1=['N of the group = ',num2str(length(grouptobeplotted))];
set(handles.text_CurrentGroupInfo,'string',info1)
info1a=['Equivalent linear distance = ',num2str(handles.clusterdistance), '  km'];
set(handles.text_LinearDistance,'string',info1a)
non_clustered= length(find(handles.datagroupindex == 0));
info4a=['Remaining non-Grouped observations: ', num2str(non_clustered)];
set(handles.text_NumberOfNonClusterData,'string',info4a)

guidata(hObject,handles)









% --- Executes on button press in pushbuttonFinishGroupSearch.
function pushbuttonFinishGroupSearch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFinishGroupSearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


delete(handles.histaux_txt)
delete(handles.histaux_lines)
delete(handles.handlehistaux)
delete(handles.histlinaux_txt)
delete(handles.histlinaux_lines)
delete(handles.handlehistlinaux)
delete(handles.htest)


            
handles.separation=handles.separation1;
if handles.groupnumber == 1
    handles.stdofgroup=handles.stdofref;
    handles.LinearDistance=handles.clusterdistance;
else
    handles.stdofgroup=[handles.stdofgroup; handles.stdofref];
    handles.LinearDistance=[handles.LinearDistance; handles.clusterdistance];
end

[ma, na]=find(handles.datagroupindex == 0);

if isempty(ma) == 0
    axes(handles.axes1)
    grouphandle=plotm(handles.inputdata(ma,1),handles.inputdata(ma,2),'co','MarkerSize',15);
    drawnow
    handles.datagroupindex(ma,1)=NaN;
end

set([handles.pushbutton_ChMembership, handles.pushbutton_NextGroup, handles.pushbuttonFinishGroupSearch],'Enable','off')

guidata(hObject,handles)





function datatemp=CreateAxes(inputdata, htoplot)
global hmap
% Test for duplicate entries on the data, and if found, then delete them

datatemp=unique(inputdata,'rows');

if length(datatemp) ~= length(inputdata)
    warning('Duplicate locations were found and removed')
    inputdata=datatemp;
end

latlim=[min(inputdata(:,1)), max(inputdata(:,1))];
lonlim=[min(inputdata(:,2)), max(inputdata(:,2))];
map_lat_extension=latlim(2)-latlim(1);
map_lon_extension=lonlim(2)-lonlim(1);
if map_lat_extension >= 150
    map_lat_lim=[-90, 90];
else
    map_lat_lim=[latlim(1)-0.1*map_lat_extension, latlim(2)+0.1*map_lat_extension];
end
if map_lon_extension >= 320
    map_lon_lim=[0 360];
else
    map_lon_lim=[lonlim(1)-0.1*map_lon_extension, lonlim(2)+0.1*map_lon_extension];
end

axes(htoplot);
if map_lat_extension >= 150 | map_lon_extension >= 320
    hmap=axesm('mollweid','MapLatLim',map_lat_lim,'MapLonLim',map_lon_lim,'MeridianLabel','on','ParallelLabel','on','Frame','on');
    setm(hmap,'MLabelLocation',30,'MLabelParallel','equator','PLabelLocation',[-80 -60 -40 -20 20 40 60 80],'PLabelMeridian','west','MLabelRound',0,'PLabelRound',0,'Grid','on')
else
    hmap=axesm('mercator','MapLatLim',map_lat_lim,'MapLonLim',map_lon_lim,'MeridianLabel','on','ParallelLabel','on','Frame','on');
    setm(hmap,'MLabelLocation',map_lon_lim,'MLabelParallel','south','PLabelLocation',map_lat_lim,'PLabelMeridian','west','MLabelRound',-2,'PLabelRound',-2)
end
tightmap



function hhist0=mkhistogram(inputdata,htoplot)

axes(htoplot);
NND=nndistance(inputdata);

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

            
function histlinear0=mklinhistogram(inputdata,htoplot)

axes(htoplot);
NND=nndistance(inputdata);

counts=hist(log(NND));
            histlinear0=histogram(NND);
            hold on
            xstd=[mean(NND,'omitnan')
            mean(NND,'omitnan') + std(NND,'omitnan')
            mean(NND,'omitnan') + 1.5*std(NND,'omitnan')
            mean(NND,'omitnan') + 1.75*std(NND,'omitnan')
            mean(NND,'omitnan') + 2*std(NND,'omitnan')];
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
            
            xlabel('NND (km)')


% --- Executes during object creation, after setting all properties.
function text3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2




% --- Executes during object creation, after setting all properties.
function popupmenuPlanet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuPlanet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text_Ellipsoid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_Ellipsoid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit_StdofRef_Callback(hObject, eventdata, handles)
% hObject    handle to edit_StdofRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_StdofRef as text
%        str2double(get(hObject,'String')) returns contents of edit_StdofRef as a double


% --- Executes during object creation, after setting all properties.
function edit_StdofRef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_StdofRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSaveFileNameNumber_Callback(hObject, eventdata, handles)
% hObject    handle to editSaveFileNameNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSaveFileNameNumber as text
%        str2double(get(hObject,'String')) returns contents of editSaveFileNameNumber as a double


% --- Executes during object creation, after setting all properties.
function editSaveFileNameNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSaveFileNameNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonSaveResults.
function pushbuttonSaveResults_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

App_number=get(handles.editSaveFileNameNumber, 'String');

longfilename=handles.filename;
shflnm=length(longfilename);

inputdata=handles.inputdata;
datagroupindex=handles.datagroupindex;
stdofgroup=handles.stdofgroup;
Lineardistances=handles.LinearDistance;
nofgroups=max(datagroupindex);
members=ones(nofgroups+1,1);
GCoordinates=zeros(nofgroups,2);

for i=1:nofgroups
    [ma na]=find(datagroupindex == i);
    members(i,1)=length(ma);
    GCoordinates(i,1)=mean(inputdata(ma,1));
    GCoordinates(i,2)=mean(inputdata(ma,2));
end
 
[ma na]=find(isnan(datagroupindex));
members(nofgroups+1,1)=length(ma);

Group_Number=1:nofgroups;
Group_Number=Group_Number';
Group_Lat=GCoordinates(:,1);
Group_Lon=GCoordinates(:,2);
Number_of_Vents=members(1:length(members)-1);
Group_Std=stdofgroup;
Group_Separation=Lineardistances;

T=table(Group_Number, Group_Lat, Group_Lon, Number_of_Vents, Group_Std, Group_Separation);


if strcmp(longfilename(shflnm-2:shflnm),'mat')
    if strcmp('NaN',App_number) ==1
        filenametobesaved=[handles.filename(1:shflnm-4),'_membership.mat'];
        filenametobesaved2=[handles.filename(1:shflnm-4),'_membership','.txt'];
    else
        filenametobesaved=[handles.filename(1:shflnm-4),'_membership',App_number,'.mat'];
        filenametobesaved2=[handles.filename(1:shflnm-4),'_membership',App_number,'.txt'];
    end
    save(filenametobesaved,'inputdata','datagroupindex','stdofgroup','Lineardistances','members','GCoordinates')
    writetable(T,filenametobesaved2,'delimiter','\t')
else
    if strcmp('NaN',App_number) ==1
        filenametobesaved=[handles.filename(1:shflnm-4),'_membership.txt'];
    else
        filenametobesaved=[handles.filename(1:shflnm-4),'_membership',App_number,'.txt'];
    end
    writetable(T,filenametobesaved,'delimiter','\t')
    
end

axes(handles.axes1)
for i=1:nofgroups
    t{i}=textm(GCoordinates(i,1),GCoordinates(i,2),num2str(i),'FontSize',30,'FontWeight','Bold','Color','y');
end

filenametobesaved3=[filenametobesaved(1:length(filenametobesaved)-4),'.jpg'];
fignew=figure;
newAxes=copyobj(handles.axes1,fignew);
set(newAxes,'Units','normalized','Position',[.10, .10, .85, .89])
saveas(fignew,filenametobesaved3)
delete(fignew)

axes(handles.axes1)
for i=1:nofgroups
    delete(t{i})
end

filenametobesaved4=[filenametobesaved(1:length(filenametobesaved)-4),'_no_numbers.jpg'];
fignew=figure;
newAxes=copyobj(handles.axes1,fignew);
set(newAxes,'Units','normalized','Position',[.10, .10, .85, .89])
saveas(fignew,filenametobesaved4)
delete(fignew)

set(handles.text_InfoSavedFile,'string','File(s) saved on curent directory')





% --- Executes on selection change in listbox_BaseImage.
function listbox_BaseImage_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_BaseImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_BaseImage contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_BaseImage


% --- Executes during object creation, after setting all properties.
function listbox_BaseImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_BaseImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_Georef.
function listbox_Georef_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_Georef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_Georef contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_Georef


% --- Executes during object creation, after setting all properties.
function listbox_Georef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_Georef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_BaseImage.
function popupmenu_BaseImage_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_BaseImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_BaseImage contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_BaseImage

global inputdata_cell 

val=get(hObject,'Value');
BaseImage_varname=handles.popupmenu_BaseImage.String{handles.popupmenu_BaseImage.Value};
handles.BaseImage=inputdata_cell.(BaseImage_varname);

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_BaseImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_BaseImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_Gref.
function popupmenu_Gref_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Gref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Gref contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Gref


global inputdata_cell 
 
val=get(hObject,'Value');
Gref_varname=handles.popupmenu_Gref.String{handles.popupmenu_Gref.Value};
handles.Gref=inputdata_cell.(Gref_varname);

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_Gref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Gref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_AddImage.
function pushbutton_AddImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_AddImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hmap

if isempty(handles.alldata) == 0
     delete(handles.alldata)
 end

map_lat_lim=handles.Gref.LatitudeLimits;
map_lon_lim=handles.Gref.LongitudeLimits;

axes(handles.axes1)

geoshow(handles.BaseImage,handles.Gref)
setm(hmap,'MapLatLim',map_lat_lim,'MapLonLim',map_lon_lim,'MLabelLocation',map_lon_lim,'MLabelParallel','south','PLabelLocation',map_lat_lim,'PLabelMeridian','west','MLabelRound',-2,'PLabelRound',-2)
tightmap
handles.alldata=plotm(handles.inputdata(:,1),handles.inputdata(:,2),'r^','MarkerFaceColor','r','MarkerSize',4);

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function text_NumberOfNonClusterData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_NumberOfNonClusterData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_Exit.
function pushbutton_Exit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
closereq


% --- Executes on button press in pushbutton_StartNewAnalysis.
function pushbutton_StartNewAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_StartNewAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
cla


if isempty(handles.Gref) == 0
    map_lat_lim=handles.Gref.LatitudeLimits;
    map_lon_lim=handles.Gref.LongitudeLimits;
    hmap=axesm('mercator','MapLatLim',map_lat_lim,'MapLonLim',map_lon_lim,'MeridianLabel','on','ParallelLabel','on','Frame','on');
    geoshow(handles.BaseImage,handles.Gref)
else
    datatemp=CreateAxes(handles.inputdata, handles.axes1);
end

handles.alldata=plotm(handles.inputdata(:,1),handles.inputdata(:,2),'r^','MarkerFaceColor','r','MarkerSize',4);

handles.stdofref=1;
set(handles.edit_StdofRef, 'String', num2str(handles.stdofref));

if isempty(handles.handlehistaux) == 0
    delete(handles.histaux_txt)
    delete(handles.histaux_lines)
    delete(handles.handlehistaux)
    delete(handles.histlinaux_txt)
    delete(handles.histlinaux_lines)
    delete(handles.handlehistlinaux)
    delete(handles.htest)
    delete(handles.grouphandle)
end

handles.groupnumber=0;
info4=['Groups identified so far: ', num2str(handles.groupnumber)];
set(handles.text_NumberofGroupInfo,'string',info4)
info1=['N of the group = '];
set(handles.text_CurrentGroupInfo,'string',info1)
info1a=['Equivalent linear distance = '];
set(handles.text_LinearDistance,'string',info1a)
non_clustered= length(handles.inputdata);
info4a=['Remaining non-Grouped observations: ', num2str(non_clustered)];
set(handles.text_NumberOfNonClusterData,'string',info4a)

set(handles.pushbuttonStartGrouping,'Enable','on')
guidata(hObject,handles)

function NND=nndistance(data)

% This function calculates the nearest neighbor distance within a group of
% data with lat, long coordinates. Earth radius is assumed by default, unless
% a global variable had been defined trhough workingplanet('planet').
% Output is a vector whose ith element indicates the distance from 
% the ith point to the datapoint closest to it.  


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

[m,n]=size(data);
separation=zeros(m-1,m);

for i=1:m
    dataaux=data;
    dataaux(i,:)=[];
    separation(:,i)=distance(data(i,1),data(i,2),dataaux(:,1),dataaux(:,2),ellipsoid);
end
NND=min(separation);


function [latpdf, lonpdf, F, Flevel]=Gausspdffinalmap(latitude, longitude, Cn, plotoption)


% Input variables are limited to:
% latitude and longitude which are the coordinates of the data, Cn is the scale
% parameter, and plotoption, which can be 0 (no plot), 1 (plot with no data) or 2 (plot with data).
% The only plot
% produced here is a default plane, with no geographic references.
% Output is formed by latpdf and longpdf which are the coordinates of the observation
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


function [G Glevel]=gausstest2(latitude,longitude,Cn,latpdf, longpdf, T, P);

% instructions used to calculate times of computation associated with the
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
    
    
    
for i=1:length(T),
    separation=distance([T(i) P(i)],[latitude, longitude],ellipsoid);
    gauss(i)=(1/(2*pi*Cn^2))*sum(exp((-1/(2*Cn^2))*separation.^2));
end


G=reshape(gauss, mT,nT);

Gmax=max(max(G));
    
Glevel=[0 .10*Gmax .20*Gmax .3*Gmax .4*Gmax .5*Gmax .6*Gmax .7*Gmax .8*Gmax .9*Gmax];
