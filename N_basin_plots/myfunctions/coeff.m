%%%functions for evaluating parallel and antiparallel (completely
%%%symmetric) modes
%%%
%parallel mode: finds all nontrivial solutions
function [condition, coefficients]=coeff(x,b,s_0,s_i_list, dimensions, n_value, tau_interval_top, tau_interval_bottom)
dim=dimensions;
w=x(end);

if (~exist("tau_interval_top", 'var'))
     % if tau_interval_top/tau_interval_bottom not specified, we default to
     % N-basin geometry meeting at one crunode
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

%mixing vectors M_{j}^{i} and vector m
M=cell(n_value-1,n_value);
m_vector=cell(n_value-1,n_value);

for l=1:n_value
    for j=1:n_value-1
        M{j,l}=zeros(2*dim, 2*dim+1);
        m_vector{j,l}=zeros(2*dim, 1);
        for m=-dim+1:dim
            partition=size(tau_interval_top);
            for d=1:partition(1)
                m_vector{j,l}(m+dim)=m_vector{j,l}(m+dim)+exp(1j*(-(m-j/n_value))*tau_interval_top{d}(l))/(1j*(-(m-j/n_value)))-exp(1j*(-(m-j/n_value))*tau_interval_bottom{d}(l))/(1j*(-(m-j/n_value)));
            end
            for n=-dim:dim
                for d=1:partition(1)
                    M{j,l}(m+dim,n+dim+1)=M{j,l}(m+dim,n+dim+1)+exp(1j*(n-(m-j/n_value))*tau_interval_top{d}(l))/(1j*(n-(m-j/n_value)))-exp(1j*(n-(m-j/n_value))*tau_interval_bottom{d}(l))/(1j*(n-(m-j/n_value))); 
                end
            end
        end
        M{j,l}(:,dim+1)=[];
    end
end

A_0_j=cell(n_value,1); %matrix A^{(0)}_{Nm-j} 
A_i=cell(n_value,1); %matrix A^{(i)} for each inner basin
for j=1:n_value
    A_0_j{j}=zeros(2*dim);
    A_i{j}=zeros(2*dim+1);
    for i=-dim+1:dim
        D_i=(-b^2+4*(b*((n_value*i-j)/(n_value*w))-((n_value*i-j)^2/n_value^2)))^(1/2);
        A_0_j{j}(i+dim,i+dim)=(-b/2+D_i*(-tan(D_i*s_0/2)^(-1))/2);
    end

    for k=-dim:dim
        D_i=(-b^2+4*(b*(k/w)-(k^2)))^(1/2);
        alpha_i=(D_i*tan(D_i*s_i_list(j)/2)/2+(b/2+abs(k)))./(D_i/2-(b/2+abs(k))*tan(D_i*s_i_list(j)/2));
        A_i{j}(k+dim+1,k+dim+1)=(-b/2+D_i*alpha_i/2);
    end
    A_i{j}(:,dim+1)=[];
    A_i{j}(dim+1,:)=[];
end

A_0_0=zeros(2*dim); %matrix A^{(0)}_{0}
for i=-dim:dim
    D_i=(-b^2+4*(b*(i/w)-(i^2)))^(1/2);
    A_0_0(i+dim+1,i+dim+1)=(-b/2+D_i*(-tan(D_i*s_0/2)^(-1))/2);
end
A_0_0(:,dim+1)=[];
A_0_0(dim+1,:)=[];

A_mat=cell(n_value-1, n_value-1);
A_hat_mat=cell(n_value-1, n_value-1);
m_hat=cell(n_value-1, n_value-1);

for j=1:n_value-1
    for i=1:n_value-1
        A_hat_mat{j,i}=(A_0_j{j}*M{j,i}-M{j,i}*A_i{i})-(A_0_j{j}*M{j,end}-M{j,end}*A_i{end})*inv(A_0_0-A_i{end})*(A_0_0-A_i{i});
        A_mat{j,i}=(-M{j,end}*inv(A_0_0-A_i{end})*(A_0_0-A_i{i})+M{j,i});
        m_hat{j,i}=m_vector{j,i};
    end
end

%composite matrices for calculating coefficients
M_mathbb=zeros((n_value-1), n_value-1);
A_mathbb=zeros(n_value-1, 2*(n_value-1)*dim);

A_mathcal=zeros(2*(n_value-1)*dim);
M_mathcal=zeros(2*(n_value-1)*dim, n_value-1);

A_mathbb_hat=zeros(2*(n_value-1)*dim);
M_mathbb_hat=zeros(2*(n_value-1)*dim, n_value-1);
for l=1:n_value-1
    for i=1:n_value-1
        sum_A_mathbb_hat=zeros(1,2*dim);
        sum_m_mathbb_hat=zeros(1, 1);
        for j=1:n_value-1
            sum_A_mathbb_hat=sum_A_mathbb_hat+m_hat{j,i}'*A_mat{j,l}/(2*pi*n_value);
            sum_m_mathbb_hat=sum_m_mathbb_hat+m_hat{j,i}'*(m_hat{j,l}-m_vector{j,end})/(2*pi*n_value); 

            A_mathbb_hat(1+2*(l-1)*dim:2*l*dim,1+2*(i-1)*dim:2*i*dim)=A_hat_mat{l,i};
            M_mathbb_hat(1+2*(l-1)*dim:2*l*dim,i)=-A_0_j{l}*(m_hat{l,i}-m_vector{l,end}); 
        end

        A_mathcal(1+2*(l-1)*dim:2*l*dim,1+2*(i-1)*dim:2*i*dim)=A_mat{l,i};
        M_mathcal(1+2*(l-1)*dim:2*l*dim,i)=(m_hat{l,i}-m_vector{l,end}); 

        A_mathbb(i,1+2*(l-1)*dim:2*l*dim)=sum_A_mathbb_hat;
        M_mathbb(i,l)=sum_m_mathbb_hat;
    end
end

%get the condition det=0
x_test=A_mathbb_hat\M_mathbb_hat;
Mat_total=A_mathbb*x_test+M_mathbb-2*pi*eye(n_value-1);
condition=det(Mat_total); 

%%%calculate coefficients from matrix by finding the eigenvectors a_0_i
a_0_i=zeros(n_value,1);
[U, S, V] = svd(Mat_total);

a_0_i(1:n_value-1)=V(:,end);
a_0_i(end)=-sum(a_0_i);

a_vector_i=zeros(2*n_value*dim,1);
a_vector_i(1:2*(n_value-1)*dim)=inv(A_mathbb_hat)*M_mathbb_hat*a_0_i(1:end-1);

a_vector_Nm=zeros(2*dim,1);
for i=1:n_value-1
   a_vector_i(1+2*(n_value-1)*dim:2*n_value*dim)=a_vector_i(1+2*(n_value-1)*dim:2*n_value*dim)-inv(A_0_0-A_i{end})*(A_0_0-A_i{i})*a_vector_i(1+2*(i-1)*dim:2*i*dim);
end
for i=1:n_value
    a_vector_Nm=a_vector_Nm+a_vector_i(1+2*(i-1)*dim:2*i*dim)/n_value;
end

%get all the coefficients of the outer basin
a_vector_0=zeros(2*n_value*dim,1);
for i=1:2*dim
    if i<=dim
        a_vector_0(n_value*i-1)=a_vector_Nm(i);
    else
        a_vector_0(n_value*i)=a_vector_Nm(i);
    end
end
a_vector_0_frac=( A_mathcal*a_vector_i(1:2*(n_value-1)*dim)+M_mathcal*a_0_i(1:end-1) )/(2*n_value*pi);
for j=1:n_value-1
    for i=1:2*dim
        term1=a_vector_0_frac(1+2*(j-1)*dim:2*j*dim);
        if i<=dim
            a_vector_0(n_value*i-j+1)=term1(i);
        else
            a_vector_0(n_value*i-j)=term1(i);
        end
    end
end

%coefficients of the vector given eigenvalue omega
coefficients={a_vector_0, a_vector_i, a_0_i};
end
