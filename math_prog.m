function varargout = math_prog(varargin)
% MATH_PROG MATLAB code for math_prog.fig
%      MATH_PROG, by itself, creates a new MATH_PROG or raises the existing
%      singleton*.
%
%      H = MATH_PROG returns the handle to a new MATH_PROG or the handle to
%      the existing singleton*.
%
%      MATH_PROG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATH_PROG.M with the given input arguments.
%
%      MATH_PROG('Property','Value',...) creates a new MATH_PROG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before math_prog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to math_prog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help math_prog

% Last Modified by GUIDE v2.5 12-May-2022 11:15:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @math_prog_OpeningFcn, ...
                   'gui_OutputFcn',  @math_prog_OutputFcn, ...
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



% --- Executes just before math_prog is made visible.
function math_prog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to math_prog (see VARARGIN)

% Choose default command line output for math_prog
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes math_prog wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = math_prog_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
n=str2double(get(handles.num_1,'string'))  %number of pionts
m=str2double(get(handles.table_1,'data'))  %matrix of two columns
t=m(:,1)    %the first column
for i=1:1:n
    
    x(i)=t(i)
end
x                      %x is avector consists of x_values from 1 to n

k=m(:,2)   %the second column
for i=1:1:n
    
    y(i)=k(i)
end
y              %y is avector consists of x_values from 1 to n


f=get(handles.list,'value')

switch f
    case 2
        
  % sum are used for storing sum of x , x^2 , y and xy to be used in solution     
 sum_X = 0;
sum_X2 = 0;
sum_Y = 0;
sum_XY = 0;

% sums are evaluated using sum function that add elements of vector
sum_X = sum(x);
sum_Y = sum(y);
sum_X2 = sum(x.^2);
sum_XY = sum(x.*y);

% soln matrix is formed then applying rref (reduced row echelon form) and
% taking it's 3rd column element which are ao and a1
soln = [n sum_X sum_Y ; sum_X sum_X2 sum_XY];
soln = rref(soln);
ao = soln(1,3);
a1 = soln(2,3);
% Evaluation a and b of the original equation
b = ao ;
a = a1 ;

% Evaluating y after linearization

 y_new= a*x + b;

avg_Y = mean(y);               % Average of y 
St = sum((y-avg_Y).^2);        % sum of difference between avg_y and y
Sr = sum((y-ao-a1*x).^2);      % sum of squared difference
r2=(St-Sr)/St ;                % coefficient of determination
r=sqrt(r2);                    % coorelation coefficient


    case 3
        
        %%Exponetial_model
    % sum are used for storing sum of x , x^2 , y and xy to be used in solution
sum_X = 0;
sum_X2 = 0;
sum_Y = 0;
sum_XY = 0;
    %since y = ax^b is linearized into ln(y) = ln(a) + bx then the
    % stored values in x and y arenot the input ones 
 
 X =x;
 Y =log(y);
 
 % sums are evaluated using sum function that add elements of vector
sum_X = sum(X);
sum_Y = sum(Y);
sum_X2 = sum(X.^2);
sum_XY = sum(X.*Y);

soln = [n sum_X sum_Y ; sum_X sum_X2 sum_XY];
soln = rref(soln);
ao = soln(1,3);
a1 = soln(2,3);

% Evaluation a and b of the original equation
a = exp(ao) ;
b = a1 ;

% Evaluating y after linearization
y_new = a*exp(x.*b);

avg_Y = mean(Y);               % Average of y 
St = sum((Y-avg_Y).^2);        % sum of difference between avg_y and y
Sr = sum((Y-ao-a1*X).^2);      % sum of squared difference
r2=(St-Sr)/St ;                % coefficient of determination
r=sqrt(r2);                    % coorelation coefficient

        
        
       
        
        
    case 4
        %%power_model
            
         % sum are used for storing sum of x , x^2 , y and xy to be used in solution
sum_X = 0;
sum_X2 = 0;
sum_Y = 0;
sum_XY = 0;
    %since y = ax^b is linearized into log(y) = log(a) + b log(x) then the
    % stored values in x and y arenot the input ones but their log
 X =log10(x);
 Y =log10(y);
 
 % sums are evaluated using sum function that add elements of vector
sum_X = sum(X);
sum_Y = sum(Y);
sum_X2 = sum(X.^2);
sum_XY = sum(X.*Y);

soln = [n sum_X sum_Y ; sum_X sum_X2 sum_XY];
soln = rref(soln);
ao = soln(1,3);
a1 = soln(2,3);

% Evaluation a and b of the original equation
a = 10^ao ;
b =  a1 ;

% Evaluating y after linearization
y_new = a*x.^b;

avg_Y = mean(Y);               % Average of y 
St = sum((Y-avg_Y).^2);        % sum of difference between avg_y and y
Sr = sum((Y-ao-a1*X).^2);      % sum of squared difference
r2=(St-Sr)/St ;                % coefficient of determination
r=sqrt(r2);                    % coorelation coefficient



        
    case 5
             %%Growthrate_model
         % sum are used for storing sum of x , x^2 , y and xy to be used in solution
sum_X = 0;
sum_X2 = 0;
sum_Y = 0;
sum_XY = 0;
    % since y = ax^b is linearized into 1/y = 1/a + b/a *(1/x) then the
    % stored values in x and y arenot the input ones but their log
 X =1./x;
 Y =1./y;
 
 % sums are evaluated using sum function that add elements of vector
sum_X = sum(X);
sum_Y = sum(Y);
sum_X2 = sum(X.^2);
sum_XY = sum(X.*Y);

soln = [n sum_X sum_Y ; sum_X sum_X2 sum_XY];
soln = rref(soln);
ao = soln(1,3);
a1 = soln(2,3);

% Evaluation a and b of the original equation
a = 1/ao ;
b =  a1/ao ;


% Evaluating y after linearization
y_new = a * x ./( b + x);

avg_Y = mean(Y);               % Average of y 
St = sum((Y-avg_Y).^2);        % sum of difference between avg_y and y
Sr = sum((Y-ao-a1*X).^2);      % sum of squared difference
r2=(St-Sr)/St ;                % coefficient of determination
r=sqrt(r2);                    % coorelation coefficient
        
end

set(handles.num_2,'string',r); % To display the correlation coeficient 
axes(handles.axes1)
plot(x,y_new);  %To plot rhe choosen model 
hold on;
plot(x,y,'.');   %To plot the pionts
title('Data Point VS Least Square Fit');
xlabel('x');
ylabel('y');




%axes(handles.axes1)
%plot(x,y)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in list.
function list_Callback(hObject, eventdata, handles)
set(handles.grid,'enable','on');
% hObject    handle to list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list


% --- Executes during object creation, after setting all properties.
function list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in grid.
function grid_Callback(hObject, eventdata, handles)
a=get(handles.grid,'value')
a
if a==1
    grid on

end
% hObject    handle to grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
set(handles.grid,'enable','off');
clc
clear 
cla
grid off
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
close
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function num_1_Callback(hObject, eventdata, handles)
% hObject    handle to num_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_1 as text
%        str2double(get(hObject,'String')) returns contents of num_1 as a double


% --- Executes during object creation, after setting all properties.
function num_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_2_Callback(hObject, eventdata, handles)
% hObject    handle to num_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_2 as text
%        str2double(get(hObject,'String')) returns contents of num_2 as a double


% --- Executes during object creation, after setting all properties.
function num_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
