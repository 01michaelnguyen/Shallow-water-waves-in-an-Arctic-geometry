
%plot parallel modes (streamfunction, velocity profile, vorticity)
function h=plot_modes(x,y,X,Y,T,S,t,omega_anti,b,s_0, n_mode, n_value, s_i_list, tau_interval_top, tau_interval_bottom, dim,basin_type,number_plots,plottype,plot_modes, plotwhat,w_solve_values,w_solve_values_anti, filter)
    %%change tau according to the basin type -- n-basin at one crunode no need
    %%to change, otherwise unwrap
    %add a sector_index for the n-basin one crunode type
    [theta, ~] = cart2pol(X, Y);  % theta is in radians, from -pi to pi
    sector_index = floor(mod(theta + pi/(n_value), 2*pi) / (2*pi/n_value)) + 1;
    if basin_type==2
        T=unwrap(T,[],2);
        if n_value==3
            T(Y>0)=T(Y>0)-2*pi;
            T=T+3*pi;
        elseif n_value==4
            T=T+3*pi;
        end
    end
    %%
    %now mark each section for each basin
    S1=S(:,:);
    S1(S<=0)=NaN;
    S1(S>s_0)=NaN;
    
    S_i_sloping=cell(n_value);
    S_i_flat=cell(n_value);
    for i=1:n_value
        S_i_sloping{i}=S(:,:);
        S_i_sloping{i}(S<s_i_list(i))=NaN;
        S_i_sloping{i}(S>0)=NaN;
    
        S_i_flat{i}=S(:,:);
        S_i_flat{i}(S>s_i_list(i))=NaN;
    end
    
    %Draw our contour!
    figure; %figure
    %writing the legend
    titleEntries = {};  % to add title descriptions
    tiledlayout(floor(number_plots^0.5),floor(number_plots^0.5), "TileSpacing","tight");
    for i=1:number_plots
    
        if plotwhat==0 
            omega=w_solve_values(1); %fix the value of omega to just the first value in the list
            omega_anti=w_solve_values_anti(1);
    
            if plot_modes==0
                time_period=linspace(0, 2*pi/omega, number_plots); %time over one period
            elseif plot_modes==1
                 time_period=linspace(0, 2*pi/omega_anti, number_plots); %time over one period
            elseif plot_modes==2
                 time_period=linspace(0, 2*pi/abs(omega-omega_anti), number_plots); %time over one period of the beating mode
            end
            t=time_period(i); %plot a certain over a given period
    
        elseif plotwhat==1
            t=0; %freeze time in one snapshot
    
            %change omega for each iteration (descending order of omega)
            if i >=length(w_solve_values)
                omega=w_solve_values(end);
            else
                omega=w_solve_values(i);
            end 
    
            %Now same but for antiparallel modes
            if i >=length(w_solve_values_anti)
                omega_anti=w_solve_values_anti(end);
            else
                omega_anti=w_solve_values_anti(i);
            end
        end
    
        %now find the coefficients 
        if basin_type==1
            [~,sol_coefficients]=coeff(omega, b, s_0, s_i_list, dim, n_value); %get parallel modes
        else
            [~,sol_coefficients]=coeff(omega, b, s_0, s_i_list, dim, n_value, tau_interval_top, tau_interval_bottom); %get parallel modes
        end
        a_vector_0=sol_coefficients{1,1};
        a_vector_i=sol_coefficients{1,2};
        a_0_i=sol_coefficients{1,3};
    
        nexttile
        %%%make outer basin s_0
        Z1=zeros(size(X));
        Z1_sectors=cell(n_value,1);
        for n=-n_value*dim:n_value*dim
            if n==0
                Z1=Z1+0;
            else
                D_n=(-b^2+4*(b*(n/(n_value*omega))-n^2/n_value^2))^(1/2);
                alpha_n_0=-tan(D_n*s_0/2)^(-1);
    
                if basin_type==1
                    for k=1:n_value
                        Z1_sectors{k}=zeros(size(X));
                        if n<0
                            Z1_sectors{k}=Z1_sectors{k}+(1-filter*abs(n)/(n_value*dim))*real(a_vector_0(n+n_value*dim+1).*exp(-b.*S1/2).*(cos(D_n.*S1/2)+alpha_n_0.*sin(D_n.*S1/2)).*exp(1j.*((T+2*(k-1)*pi)*n/n_value-omega*t)));
                        else
                            Z1_sectors{k}=Z1_sectors{k}+(1-filter*abs(n)/(n_value*dim))*real(a_vector_0(n_value*dim+n).*exp(-b.*S1/2).*(cos(D_n.*S1/2)+alpha_n_0.*sin(D_n.*S1/2)).*exp(1j.*((T+2*(k-1)*pi)*n/n_value-omega*t)));
                        end
                        Z1_sectors{k}=Z1_sectors{k}.*(sector_index==k);
                        Z1_sectors{k}(isnan(Z1_sectors{k}))=0;
                        Z1=Z1+Z1_sectors{k};
                    end
                else
                    if n<0
                        Z1=Z1+(1-filter*abs(n)/(n_value*dim))*real(a_vector_0(n_value*dim+n+1).*exp(-b.*S1/2).*(cos(D_n.*S1/2)+alpha_n_0.*sin(D_n.*S1/2)).*exp(1j.*(T*n/n_value-omega*t)));
                    else
                        Z1=Z1+(1-filter*abs(n)/(n_value*dim))*real(a_vector_0(n_value*dim+n).*exp(-b.*S1/2).*(cos(D_n.*S1/2)+alpha_n_0.*sin(D_n.*S1/2)).*exp(1j.*(T*n/n_value-omega*t)));
                    end
                end
            end
        end
        Z1(isnan(Z1))=0;
    
        %%
        %%%add up each inner basin individually from i=1 to i=N
        Z_i_sloping=cell(n_value);
        Z_i_flat=cell(n_value);
    
        for j=1:n_value
            Z_i_sloping{j}=zeros(size(X));
            Z_i_flat{j}=zeros(size(X));
            for n=-dim:dim
                if n==0
                    Z_i_sloping{j}=Z_i_sloping{j}+real(a_0_i(j)*exp(-1j.*omega*t));
                    Z_i_flat{j}=Z_i_flat{j}+real(a_0_i(j)*exp(-1j*omega*t));
                else
                    vector_basin=a_vector_i(1+2*(j-1)*dim:2*j*dim);
    
                    D_n=(-b^2+4*(b*n/(omega)-(n^2)))^(0.5);
                    h_n=((D_n/2)*tan(D_n*s_i_list(j)/2)+(b/2+abs(n)))/(D_n/2-(b/2+abs(n))*tan(D_n*s_i_list(j)/2));
                    b_vector_i=vector_basin.*exp((-b/2-abs(n))*s_i_list(j))*(h_n*sin(D_n*s_i_list(j)/2)+cos(D_n*s_i_list(j)/2));
                    if n<0
                        Z_i_sloping{j}=Z_i_sloping{j}+(1-filter*abs(n)/dim)*real(vector_basin(dim+n+1).*exp(-b.*S_i_sloping{j}/2).*(cos(D_n.*S_i_sloping{j}/2)+h_n.*sin(D_n.*S_i_sloping{j}/2)).*exp(1j.*(T*n-omega*t)) );
                        Z_i_flat{j}=Z_i_flat{j}+(1-filter*abs(n)/dim)*real((b_vector_i(dim+n+1).*exp(1j*n*T)).*exp(abs(n).*S_i_flat{j})*exp(-1j.*omega.*t));
                    else
                        Z_i_sloping{j}=Z_i_sloping{j}+(1-filter*abs(n)/dim)*real(vector_basin(dim+n).*exp(-b.*S_i_sloping{j}/2).*(cos(D_n.*S_i_sloping{j}/2)+h_n.*sin(D_n.*S_i_sloping{j}/2)).*exp(1j.*(T*n-omega*t)) );
                        Z_i_flat{j}=Z_i_flat{j}+(1-filter*abs(n)/dim)*real((b_vector_i(dim+n).*exp(1j*n*T)).*exp(abs(n).*S_i_flat{j})*exp(-1j.*omega.*t));
                    end
                end
            end
            Z_i_sloping{j}(isnan(Z_i_sloping{j}))=0;
            Z_i_flat{j}(isnan(Z_i_flat{j}))=0;
        end
        
        %%now add the symmetrical solutions and repeat for the outer basin
        [Z1_anti,Z2_anti,Z3_anti]=plot_antiparallel(X,T,S,t,omega_anti,b,s_0, n_mode, n_value, s_i_list);
        Z_anti=Z1_anti+Z2_anti+Z3_anti;
    
        %calculate the streamfunction and velocities of each basin then add
        %them up
        u=zeros(size(X)); %velocity u
        v=zeros(size(X)); %velocity v
        T_j=cell(n_value,1); %tau section corresponding to basin j
    
        dx = x(2)-x(1);
        dy = y(2)-y(1);
        Z=Z1;
        if basin_type==1
            for j=1:n_value
                T_j{j}=(sector_index==j);
                Z=Z+(Z_i_sloping{j}+Z_i_flat{j}).*(sector_index==j);
                
                %plot velocities 
                [z1i_sloping_x,z1i_sloping_y]=gradient(Z1+Z_i_sloping{j}, dx, dy); %get gradients of inner basin streamfunction
                [zi_flat_x,zi_flat_y]=gradient(Z_i_flat{j}, dx, dy);
        
                u_i=-z1i_sloping_y./(exp(-b.*(S-s_i_list(j))))-zi_flat_y./(exp(-b*(s_i_list(1)-s_i_list(j))));
                u_i=u_i.*(sector_index==j);
                v_i=z1i_sloping_x./(exp(-b.*(S-s_i_list(j))))+zi_flat_x./(exp(-b*(s_i_list(1)-s_i_list(j))));
                v_i=v_i.*(sector_index==j);
                u=u+u_i;
                v=v+v_i;
            end
        else
            for j=1:n_value
                T_j{j}=zeros(size(T));
                for d=1:partition
                    T_j{j}=T_j{j}+(T<tau_interval_top{d}(j) & T>tau_interval_bottom{d}(j));
                    Z=Z+(Z_i_sloping{j}+Z_i_flat{j}).*(T<tau_interval_top{d}(j) & T>tau_interval_bottom{d}(j));
                    
                    %plot velocities 
                    [z1i_sloping_x,z1i_sloping_y]=gradient(Z1+Z_i_sloping{j}, dx, dy); %get gradients of inner basin streamfunction
                    [zi_flat_x,zi_flat_y]=gradient(Z_i_flat{j}, dx, dy);
            
                    u_i=-z1i_sloping_y./(exp(-b.*(S-s_i_list(j))))-zi_flat_y./(exp(-b*(s_1-s_2)));
                    u_i=u_i.*(T<tau_interval_top{d}(j) & T>tau_interval_bottom{d}(j));
                    v_i=z1i_sloping_x./(exp(-b.*(S-s_i_list(j))))+zi_flat_x./(exp(-b*(s_1-s_2)));
                    v_i=v_i.*(T<tau_interval_top{d}(j) & T>tau_interval_bottom{d}(j));
                    u=u+u_i;
                    v=v+v_i;
                end
            end
        end
    
        %work out velocities taking in account height => u . h = z x grad(psi), h
        %is our height depending on the value of sigma (which this takes int account)
    
        %velocities of the parallel mode
        %use this kinetic energy (parallel mode) to normalise energy of each mode to 1
        KE = 0.5 * sum(sum(u.^2 + v.^2)) * dx * dy; %compute the kinetic energy over the grid
        
        [z12_x_anti,z12_y_anti]=gradient(Z1_anti+Z2_anti, dx, dy); %get gradients of inner basin streamfunction
        [z3_x_anti,z3_y_anti]=gradient(Z3_anti, dx, dy);
        u_anti=-z12_y_anti./(exp(-b.*(S-s_i_list(2))))-z3_y_anti./1;
        v_anti=z12_x_anti./(exp(-b.*(S-s_i_list(2))))+z3_x_anti./1;
    
        %use this kinetic energy (antiparallel mode) to normalise energy of each mode to 1
        KE_anti = 0.5 * sum(sum(u_anti.^2 + v_anti.^2)) * dx * dy; %compute the kinetic energy over the grid
    
        if plottype==0 && plot_modes==0
            h=contour(X, Y, Z, 50);
            clim([min(Z,[], "all"), max(Z, [], "all")])
        elseif plottype==0 && plot_modes==1
            h=contour(X, Y, Z_anti, 50);
            clim([min(Z_anti,[], "all"), max(Z_anti, [], "all")])
        elseif plottype==0 && plot_modes==2
            h=contour(X, Y, Z/KE^0.5+Z_anti/KE_anti^0.5, 50);
            clim([min(Z/KE^0.5+Z_anti/KE_anti^0.5,[], "all"), max(Z/KE^0.5+Z_anti/KE_anti^0.5, [], "all")])
        elseif plottype==2 && plot_modes==0
            vort=divergence(X,Y,u/KE^0.5,v/KE^0.5); %plot vorticity if parallel
            h=contour(X,Y,vort,50);
        elseif plottype==2 && plot_modes==1
            vort=divergence(X,Y,u_anti/KE_anti^0.5,v_anti/KE_anti^0.5); %plot vorticity if antiparallel
            h=contour(X,Y,vort,50);
        elseif plottype==2 && plot_modes==2
            vort=divergence(X,Y,u/KE^0.5+u_anti/KE_anti^0.5,v/KE^0.5+v_anti/KE_anti^0.5);     %plot vorticity if transferring
            h=contour(X,Y,vort,50);
        elseif plottype==1 && plot_modes==0
            h=quiver(X, Y, u/KE^0.5, v/KE^0.5); %%plot velocity field if parallel
        elseif plottype==1 && plot_modes==1
            h=quiver(X, Y, u_anti/KE_anti^0.5, v_anti/KE_anti^0.5); %%plot velocity field if antiparallel
        elseif plottype==1 && plot_modes==2
            h=quiver(X, Y, u/KE^0.5+u_anti/KE_anti^0.5, v/KE^0.5+v_anti/KE_anti^0.5); %%plot velocity field if transferring
        end
    
        hold on
        for j=1:n_value
            contour(X,Y, S.*T_j{j}, [s_i_list(j), s_i_list(j)], "k-.") %get contours of the basin itself
        end
        contour(X,Y, S, [0, 0], "k:") %get contours of the basin itself
        contour(X,Y, S,[s_0, s_0], "k-.")
        axis off
        hold off
    
        %make diverging colormap
        cmocean('balance','pivot',0)
        %title(titleEntries{end}, "interpreter", "latex", "fontsize", 14); %plot title
    end
    

end