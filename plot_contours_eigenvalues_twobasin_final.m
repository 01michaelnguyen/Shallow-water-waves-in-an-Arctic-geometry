%Plot the free mode eigenvalue solutions of the linearised rigid-lid Shallow water equations
%on the two-basin geometry
options=optimset('MaxIter',1e10,'MaxFunEvals',1e10, "TolFun", 1e-10, "TolX", 1e-20);
addpath('myfunctions')


%%%%
%%%% PLOT CONFIGURATION SETTINGS
%%%%


%%configure what are we plotting
plottype=0; %plottype=0 means plot psi, plottype=1 means plot the velocity field, plottype=2 means plot vorticity
plot_modes=0; 
%plot_modes=0 means plot just the parallel type modes, 
% plot_modes=1 means plot the antiparallel type modes (only applicable if s_1=s_2), 
% plot_modes=2 means plot both at the same time to get the beating modes (when frequencies are close together)
plotwhat=0; 
%plotwhat=0 means plot ONE (the highest omega) mode propagates in one period 
%plotwhat=1 means plot different modes in a given time snapshot

%%n-basin problem:
n_value=2; %set number of basins (n_value=2 is the 2 basin problem)
basin_type=1; 
%basin_type=1 means joined basins meeting at one point
%basin_type=2 means basins aligned along a line (only works for n_value=3,4 basins)

%%change resolution of plots
res=200; %resolution
filter=0; %filter=0 means no filter, filter=1 is the fejer filter, anything inbetween is varying levels of filter strength

%%change geometric parameters of the basin
b=1; %basin slope
s_0=0.6; %outer basin size
s_1=-1.5; %inner basin 1 maximum depth
s_2=-1.5; %inner basin 2 maximum depth
s_i_list=[s_1,s_2];

%%specify the azimuthal mode number of the antiparallel solution (n_mode>0)
n_mode=2; 

%%set parameters of solver and eigenvalue guesses
%%dimensions of matrix solver (higher=more accurate plots)
dim=200;
%plot the minimum and maximum (guess) value of omega for the highest mode and lowest mode in that interval 
omega_min=0.17; %0.50, 0.52, s0=3.5
omega_max=0.366;

number_plots=9; %plot how many frames per period or number of modes shown on display (should be a square number)



%%%%% 
%%%%% END OF PLOT CONFIGURATION SETTINGS
%%%%%


x=linspace(-2,2,res); %size of x-grid
y=linspace(-1.5,1.5, res); %size of y-grid
[X,Y]=meshgrid(x,y);
Z=X+1j.*Y;

if basin_type==1
    S=real(log(Z.^n_value-1));
    T=imag(log(Z.^n_value-1));
    
    partition=1;
    tau_interval_top=cell(1,1);
    tau_interval_bottom=cell(1,1);
    tau_interval_top{1}=zeros(n_value,1);
    tau_interval_bottom{1}=zeros(n_value,1);
    for i=1:n_value
        tau_interval_top{1}(i)=(2*i-1)*pi;
        tau_interval_bottom{1}(i)=(2*i-3)*pi;
    end
end
%%solve the eigenvalue problem to get omega
w_solve_values=[]; 
w_solve_values_anti=[]; 

eigen_vector_list=[];
for n=omega_max:-0.001:omega_min
    try 
        x_eigenvalue_init=n; %guess value of omega
        if basin_type==1
            fun=@(x) coeff(x, b, s_0, s_i_list, dim, n_value); %get parallel modes
        else
            fun=@(x) coeff(x, b, s_0, s_i_list, dim, n_value, tau_interval_top, tau_interval_bottom); %get parallel modes
        end
        antifun=@(w) antiparallel(w,b,s_0,s_1,n_mode); %get antiparallel modes (if applicable)

        %solve parallel mode
        w_solve=fsolve(fun, x_eigenvalue_init, options);
        %solve antiparallel mode
        w_solve_anti=fsolve(antifun, x_eigenvalue_init, options);

        if ismembertol(abs(w_solve(end)),abs(w_solve_values), 10^-6)==0 && abs(w_solve)<1
            w_solve_values=[w_solve_values, w_solve(end)];
        end
        if ismembertol(w_solve_anti(end),w_solve_values_anti, 10^-6)==0
            w_solve_values_anti=[w_solve_values_anti, w_solve_anti(end)];
        end
    catch
        disp("no solution found")
        continue
    end
    if plotwhat==0 && length(w_solve_values)==1 && length(w_solve_values_anti)==1
        break
    else
        continue
    end
end

%%list of omega values and sorted in descending order
w_solve_values=sort(w_solve_values, "descend");
w_solve_values_anti=sort(w_solve_values_anti, "descend");

%plot modes
h=plot_contours(x,y,X,Y,T,S,b,s_0, n_mode, n_value, s_i_list, ...
    dim,basin_type,number_plots,plottype,plot_modes, plotwhat,w_solve_values,w_solve_values_anti, ...
    filter,partition, tau_interval_top, tau_interval_bottom);