%plot antiparallel mode given n_mode and omega
function [Z1_anti,Z2_anti,Z3_anti]=plot_antiparallel(X,T,S,t,omega_anti,b,s_0, n_mode, n_value, s_i_list)
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

    %%now add the symmetrical solutions and repeat
    %outer basin
    Z1_anti=zeros(size(X));
    Z2_anti=zeros(size(X));
    Z3_anti=zeros(size(X));
    
    D_n=(-b^2+4*(b*n_mode/(omega_anti)-(n_mode^2)))^(0.5);
    v_n=-tan(D_n*s_0/2)^(-1);
    d_n=((D_n/2)*tan(D_n*s_i_list(1)/2)+(b/2+n_mode))/(D_n/2-(b/2+n_mode)*tan(D_n*s_i_list(1)/2));
    c_vector_anti=exp((-b/2-n_mode)*s_i_list(1))*(d_n*sin(D_n*s_i_list(1)/2)+cos(D_n*s_i_list(1)/2));

    Z1_anti=Z1_anti+real(exp(-b.*S1/2).*(cos(D_n.*S1/2)+v_n.*sin(D_n.*S1/2)).*exp(1j.*(T*n_mode-omega_anti*t)));
    Z2_anti=Z2_anti+real(exp(-b.*S_i_sloping{1}/2).*(cos(D_n.*S_i_sloping{1}/2)+d_n.*sin(D_n.*S_i_sloping{1}/2)).*exp(1j.*(T*n_mode-omega_anti*t)));
    Z3_anti=Z3_anti+real(c_vector_anti*exp(1j*n_mode*T).*exp(n_mode.*S_i_flat{1}-1j.*omega_anti.*t));
    
    Z1_anti(isnan(Z1_anti))=0;
    Z2_anti(isnan(Z2_anti))=0;
    Z3_anti(isnan(Z3_anti))=0;
end