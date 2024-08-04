%%
% Numerical Methods

% Linear regression

% Comment for GUI is inside GUI function 'math_prog.m'
math_prog;
%%
% Jacobi Method
clc; 
clear;
format LONG;
A = zeros(3,3);
B = zeros(1,3);
for i = 1 :3
    for j = 1:3
     A(i,j) = input(['Enter a' num2str(i) num2str(j) ' = ' ]);
    end
    B(i,1) = input(['Enter b' num2str(i) ' = ' ]);
end

x1 = input('Enter x1 = '); 
x2 = input('Enter x2 = '); 
x3 = input('Enter x3 = '); %in initial guess size(X) = size(B)
option = menu ('Choose which terminate Method','Error','Ilterations');

n = input('Enter Number of ilterations = '); % Maximum itration 

%Method_________________________________________________________________


sum = 0;
flag = zeros(1,3);
k =  0;
A_B = [A B];
A_B_new =  A_B;
while( k < 6)
    A = A_B_new([1 2 3],[1 2 3]);
for i = 1:3
   sum = 0;
    for j = 1 : 3
        if(i == j)
            continue;
        end
        sum = sum + abs(A(i,j));
    end
    if(abs(A(i,i)) > sum)
        flag(1,i) = 1;
    end
end
if(flag == [1 1 1] )
    break;
end
k = k+1;
switch(k)
    case 1
        A_B_new([1 2 3],:) = A_B([1 3 2],:);
    case 2
        A_B_new([1 2 3],:) = A_B([2 1 3],:);
    case 3
        A_B_new([1 2 3],:) = A_B([2 3 1],:);
    case 4
        A_B_new([1 2 3],:) = A_B([3 2 1],:);
    case 5
        A_B_new([1 2 3],:) = A_B([3 1 2],:);
end

if(k == 6)
    fprintf('Matrix cant be strictly dominant\n');
end
end
A = A_B_new(1:3,1:3);
B = A_B_new(:,4);
%Output___________________________________________________________
fprintf('ilteration\tx1\tx2\tx3\n');
fprintf('------------------------------------------------');

for i = 1:n
   x1_new = (B(1,1) - A(1,2)*x2 - A(1,3)*x3) / A(1,1);
   x2_new = (B(2,1) - A(2,1)*x1 - A(2,3)*x3) / A(2,2);
   x3_new = (B(3,1) - A(3,1)*x1 - A(3,2)*x2) / A(3,3);
   x1 = x1_new;
   x2 = x2_new;
   x3 = x3_new;
   fprintf('\n%d         \t %.5f\t %.5f\t %.5f\n',i,x1, x2,x3 );
end


%%
% Bisection Method
% Method is used to evaluate zeroes of input function bonded by limits

clc; 
clear;

format long
% this line request user to enter function but user must follow some rules:
% 1- function must be in one variable which is x
% 2- variable x must be in closed bracket
% 3- operator are the same as used in matlab code
input_function_by_user = input('Enter function in x (note use parathess with x e.g: (x) ): ','s');
% input function is converted from string to function stored in f_x
f_x = str2func(['@(x)' input_function_by_user ]);

% user is asked to enter limits of function
A = input('Enter start at a = '); %left limit
B = input('Enter end at b = ');  % Right limit

% menu is used to offer user to choose termination method either error or
% iterations returning value corresponding to the choice in option.
option = menu ('Choose which terminate Method','Error','Ilterations');
switch(option)
    case 1
         error   = input('Enter Maximum Error : '); % maximum Error
         n = ceil(log2((B-A)/error));
    case 2 
         n = input('Enter Number of ilterations = '); % Maximum itration 
         error = (B - A)/(2^n);
end 

List = zeros(n,7);
%Method_________________________________________________________________
fa = f_x(A); 
fb = f_x(B);
for i = 1:n
    r  = (A+B)/2;   
    fr = f_x(r);
    List(i,:) = [i, A, fa, B, fb, r, fr]; %for illustration
    if (i~=1 && abs(r-r0) <= error) 
        break;
    end
    if  fa*fr < 0,    B = r;  fb = fr;
    else,             A = r;  fa = fr;
    end
    r0 = r;
end

%Illustration___________________________________________________________
fprintf(['%3s',repmat('%11s',[1,6]),'\n'],...
    'itr','a','f(a)','b','f(b)','r','f(r)');
div = ['---',repmat('-----------',[1,6]), '\n'];        fprintf(div);
fprintf(['%3.0f', repmat('%11.4g',[1,6]), '\n'],List.'); fprintf(div);
fprintf('x = %10.4g\n',r);
fprintf('Error = %f\n',error);
%%
% Newton-Raphson method
% Method is used to evaluate zeroes of function input

clear;
clc;

% format function is used to change number of number after decimal
% displayed, it is chosen long to make it able to display 5 numbers.
format long;


% this line request user to enter function but user must follow some rules:
% 1- function must be in one variable which is x
% 2- variable x must be in closed bracket
% 3- operator are the same as used in matlab code
input_function_by_user = input('Enter function in x (note use parathess with x e.g: (x) ): ','s');
% input function is converted from string to function stored in f_x
f_x = str2func(['@(x)' input_function_by_user ]);

% syms fucntion is used to declare that x is varible in function that will
% be derivatived.
syms x;
% Evaluating differentiation of input function
df_x = eval(['@(x)' char(diff(f_x(x)))]);

% User is asked to enter initial x
x0 = input('Enter initial x = ');

% menu is used to offer user to choose termination method either error or
% iterations returning value corresponding to the choice in option.
option = menu ('Choose which terminate Method','Error','Iterations');
switch(option)
    case 1
         error   = input('Enter Maximum Error : '); % maximum Error
         i=0;
         fprintf('(0) X = %.5f\n',x0);
         x = x0 - (f_x(x0)/df_x(x0));

         while( abs(x-x0) > error)
             i=i+1;
             fprintf('(%d) X ',i);
             fprintf('= %.5f',x);
             fprintf('\t  E = %.5f\n',abs(x-x0));
             x0 = x;
             x = x0 - (f_x(x0)/df_x(x0));
         end
    case 2 
         n = input('Enter Number of iterations = '); % Maximum itrations
         fprintf('(0) X = %.5f\n',x0);
          for i = 1:n  
             
             x = x0 - (f_x(x0)/df_x(x0));
             if(i==n)
                 break;
             end
             fprintf('(%d) X ',i);
             fprintf('= %.5f',x);
             fprintf('\t  E = %.5f\n',abs(x-x0));
             x0 = x;  
          end
          i=i-1;
         
end 
% Displaying Output
 fprintf('(%d) X = %.5f',i+1,x);
 fprintf('\t  E = %.5f\n',abs(x-x0));
 fprintf('\nX = %.5f\n',x);
%%
% Trapezoidal method
% Method is used to calculate approxmiately value of definite integration

clear;
clc;

% User is asked to enter:
% initial point of integeration a.
% final point of integration b.
A = input('Enter start at a = ');
B = input('Enter end at b = ');

% User is asked to either enter step size or segments number.
% then the other parameter is calculated
option = menu ('Choose which Method','Step Size','Segments Number');
switch(option)
    case 1
        h = input('Enter Step size h = ');
        n=ceil((B-A)/h) ;
    case 2
        n = input('Enter number of segments n = ');
        h = (B-A)/n ;
end

% this line request user to enter function but user must follow some rules:
% 1- function must be in one variable which is x
% 2- variable x must be in closed bracket
% 3- operator are the same as used in matlab code
input_function_by_user = input('Enter function in x (note use parathess with x e.g: (x) ): ','s');
% input function is converted from string to function stored in f_x
f_x = str2func(['@(x)' input_function_by_user ]);

% Evaluating f at b and f at a
f_b = f_x(B);
f_a = f_x(A);

% f_n is vector used to store f_x at each step to be summed later
f_n = [];
for i = 1 : n-1 
    f_n = [f_n f_x(A + i*h)];
end

% Evaluating sum
sum_points = sum(f_n);

% Applying Trapeziodal Method rule
integration = (h/2)*(f_a + f_b + 2*sum_points);

% displaying the result of integration using Trapeziodal rule
fprintf(['\napproxmiate Integration = ' num2str(integration) '\n']);

%%
% Simpson's 1/3 rule
% Method is used to calculate approxmiately value of definite integration

clear;
clc;

% n is number of segment and it is initailized as 1 then it is checked if
% it is even or not as this method require odd number of points thus even
% number of segments
n=1;

% while loop check if n is even or not.
% in case of odd, repeat.
% in case of even, continue the rest of the code.
while(mod(n,2)==1)
    
% User is asked to enter:
% initial point of integeration a.
% final point of integration b.
A = input('Enter start at a = ');
B = input('Enter end at b = ');

% User is asked to either enter step size or segments number.
% then the other parameter is calculated
option = menu ('Choose which Method','Step Size','Segments Number');
switch(option)
    case 1
        h = input('Enter Step size h = ');
        n=ceil((B-A)/h) ;
    case 2
        n = input('Enter number of segments n = ');
        h = (B-A)/n ;
end


% n is evaluated as ceil of the (b-a)/h as n must be even integer

% Message is displayed in case entered h results in odd n
if(mod(n,2)==1)
    fprintf('\nNumber of segment isnot even number, Enter different h and try again \n')
end
end

% this line request user to enter function but user must follow some rules:
% 1- function must be in one variable which is x
% 2- variable x must be in closed bracket
% 3- operator are the same as used in matlab code
input_function_by_user = input('Enter function in x (note use parathess with x e.g: (x) ): ','s');
% input function is converted from string to function stored in f_x
f_x = str2func(['@(x)' input_function_by_user ]);

% Evaluating f at b and f at a
f_b = f_x(B);
f_a = f_x(A);

% this vector is used to store f_n elements for evaluating odd and even sum
f_n = [];

% Loop used for evaluating f at each step
for i = 1 : n-1 
    f_n = [f_n f_x(A + i*h)];
end

% Evaluating even and odd sums
sum_odd = sum(f_n(1:2:n-1));
sum_even =  sum(f_n(2:2:n-1));

% Applying Simpsons 1/3 rule
integration = (h/3)*(f_a + f_b + 4*sum_odd + 2*sum_even);

% displaying the result of integration using Simpson's 1/3 rule
fprintf(['\n approxmiate Integration = ' num2str(integration) '\n']);
%%
% Euler's method
% Method used to estimate solution of ODE approximately

clear;
clc;

% this line request user to enter function but user must follow some rules:
% 1- function must be in one or two variable which are x and y
% 2- variables x and y must be in closed bracket
% 3- operator are the same as used in matlab code
input_function_by_user = input('Enter function in x (note use parathess with x e.g: (x) ): ','s');
% input function is converted from string to function stored in f_x
f_x = str2func(['@(x,y)' input_function_by_user ]);

% User is asked to enter:
% step size h.
% Starting point of the curve x_start.
% ending point of the curve x_end.
h = input('Enter step size h = ');
x_start = input('Enter initial x = ');
x_end = input('Enter final x = ');

% n is evaluated as ceil of the (x_end - x_start)/h as n must be integer
n = ceil((x_end - x_start)/h );

% user is asked to enter initial y at x = 0
y_o = input(['Enter at x = ' num2str(x_start) ' , y = ']);

% y_old is initialized as y_o so that it will be used in first loop
y_old = y_o;

% y vector will store all values of y as result of Euler's method
y = [y_o];

% x is vector used as independent variables in Euler's Method
x = x_start:h:x_end;
for i = 1 : n
y_new = y_old + h* f_x(x(i),y_old);
y_old = y_new;
y = [y y_new];
end

% Displaying final result
disp(['y = ' num2str(y)]);

% Ploting the final result
figure;
plot(x,y);
%%
% Heun's method

clear;
clc;
% this line request user to enter function but user must follow some rules:
% 1- function must be in one or two variable which are x and y
% 2- variables x and y must be in closed bracket
% 3- operator are the same as used in matlab code
input_function_by_user = input('Enter function in x (note use parathess with x e.g: (x) ): ','s');
% input function is converted from string to function stored in f_x
f_x = str2func(['@(x,y)' input_function_by_user ]);

% User is asked to enter:
% step size h.
% Starting point of the curve x_start.
% ending point of the curve x_end.
h = input('Enter step size h = ');
x_start = input('Enter initial x = ');
x_end = input('Enter final x = ');

% n is evaluated as ceil of the (x_end - x_start)/h as n must be integer
n = ceil((x_end - x_start)/h );


% user is asked to enter initial y at x = 0
y_o = input(['Enter at x = ' num2str(x_start) ', y = ']);

% y_old is initialized as y_o so that it will be used in first loop
y_old = y_o;

% y vector will store all values of y as result of Heun's method
y = [y_o];

% x is vector used as independent variables in Heun's Method
x = x_start : h : x_end;
for i = 1 : n
y_star = y_old + h* f_x(x(i),y_old); % Predictor 
y_new = y_old + (h/2)*( f_x(x(i+1),y_star) + f_x(x(i) ,y_old)); % Corrector

% y_old is then become y_new in preperation for next step
y_old = y_new;

% Storing value of y_new in y for final result
y = [y y_new];
end

% Displaying final result
disp(['y = ' num2str(y)]);

% Ploting the final result
figure;
plot(x,y);