%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NAME OF FILE: POD_multi_layer_2D_ppwd_complex_epsr_var_for_thesis.m
%
% PURPOSE: POD method for parallel plate waveguide
%
% Written by Gladwin Jos
%
% Date : 09/08/2018 - lossless dielectric
% Date : 26/08/2018 - lossy dielectric
% Date : 23/09/2018 - implemented in multi layer
% Date : 14/07/2019 - reduced basis method implementation
% Date : 04/11/2019 - Extended it for parallel plate waveguide
% Date : 20/07/2020 -  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Reference: 
% "Fast Solution of High Stochastic Dimensional EM Problems Using Proper Orthogonal
% Decomposition", IEEE Microwave and Wireless components letters

clc
clear all; % Clear all variables
close all; % Close all graphics 
tic;
Y0=0;
Y1=3.5;
X0=0;
X1=25;
Y2=2.5;  
mesh_file_name='2d_parallel_plate_wgd_multilayer_dielectric_brick_25_305_height_205_adaptive_Ne_525.msh';
[node_list one_node_point two_node_line three_node_triangle ...
        four_node_tetrahedron ]=parse_gmsh_file(mesh_file_name);
n=three_node_triangle;
No_of_dielectric_layer=5; 
No_of_rv=2*No_of_dielectric_layer;
x=node_list(:,1);
y=node_list(:,2); 
Ne=length(three_node_triangle);
No_of_nodes=length(node_list);  
mu_0=4*pi*10^(-7);
epsilon_r_1_dash=2;
epsilon_r_2_dash=4;
epsilon_r_3_dash=3;
epsilon_r_4_dash=8;
epsilon_r_5_dash=6;
epsilon_r_1_dd=0.8;
epsilon_r_2_dd=0.12;
epsilon_r_3_dd=0.3;
epsilon_r_4_dd=0.22;
epsilon_r_5_dd=0.1;
freq=3.9*10^9; 
epsilon_r_air=1;
epsilon_0=8.854187817*10^(-12);
k0=2*pi*freq*sqrt(mu_0*epsilon_0);
E0=1; % amplitude


count1=1;
count2=1;
count_left=1;
count_right=1;
count_down=1;
count_top=1;
for ii=1:No_of_nodes
    if y(ii)==Y0 || y(ii)==Y1
       bc_node_PEC(count1)=ii;
       bc_value_PEC(count1)=0;
       count1=count1+1;
       if y(ii)==Y0
            bc_node_ABC_down(count_down)=ii;
            count_down=count_down+1;
       else
           bc_node_ABC_top(count_top)=ii;
           count_top=count_top+1; 
       end
    end
    if x(ii)==X0 || x(ii)==X1
       bc_node_ABC(count2)=ii;
       count2=count2+1;
       if x(ii)==X0
            bc_node_ABC_left(count_left)=ii;
            count_left=count_left+1;
       else
           bc_node_ABC_right(count_right)=ii;
           count_right=count_right+1; 
       end
    end
end

element_left_count=1;
element_right_count=1;
for e=1:Ne
    node1=n(e,1);
    node2=n(e,2);
    node3=n(e,3);
    if (x(node1)==X0 && x(node2)==X0)||(x(node1)==X0 && x(node3)==X0) ...
            ||(x(node2)==X0 && x(node3)==X0)
        element_left_bc(element_left_count,1)=e;
        element_left_count=element_left_count+1;
    elseif (x(node1)==X1 && x(node2)==X1)||(x(node1)==X1 && x(node3)==X1)...
            ||(x(node2)==X1 && x(node3)==X1)        
        element_right_bc(element_right_count,1)=e;
        element_right_count=element_right_count+1;
    end
end

% defining dielectric dimensin
% for dielectric dimension
% for layer 1
X_left_1=10;
X_right_1=X_left_1+1;
Y_bottom_1=Y0;
Y_top_1=Y2; 
% for layer 2
X_left_2=X_right_1;
X_right_2=X_left_2+1;
Y_bottom_2=Y0;
Y_top_2=Y2; 
% for layer 3
X_left_3=X_right_2;
X_right_3=X_left_3+1;
Y_bottom_3=Y0;
Y_top_3=Y2; 
% for layer 4
X_left_4=X_right_3;
X_right_4=X_left_4+1;
Y_bottom_4=Y0;
Y_top_4=Y2; 
% for layer 5
X_left_5=X_right_4;
X_right_5=X_left_5+1;
Y_bottom_5=Y0;
Y_top_5=Y2; 

% in meter
X_left_1=X_left_1/100;
X_left_2=X_left_2/100;
X_left_3=X_left_3/100;
X_left_4=X_left_4/100;
X_left_5=X_left_5/100;
X_right_1=X_right_1/100;
X_right_2=X_right_2/100;
X_right_3=X_right_3/100;
X_right_4=X_right_4/100;
X_right_5=X_right_5/100;
Y_bottom_1=Y_bottom_1/100;
Y_bottom_2=Y_bottom_2/100;
Y_bottom_3=Y_bottom_3/100;
Y_bottom_4=Y_bottom_4/100;
Y_bottom_5=Y_bottom_5/100;
Y_top_1=Y_top_1/100;
Y_top_2=Y_top_2/100;
Y_top_3=Y_top_3/100;
Y_top_4=Y_top_4/100;
Y_top_5=Y_top_5/100;
x=node_list(:,1)/100; % in meters
y=node_list(:,2)/100; % in meters
x=roundn(x,-12);
y=roundn(y,-12);
epsilon_r_final=zeros(Ne,1);
for e=1:Ne  
    node1=n(e,1);
    node2=n(e,2);
    node3=n(e,3); 
    if x(node1)>=X_left_1 && x(node1)<=X_right_1 && ...
          x(node2)>=X_left_1 && x(node2)<=X_right_1 && ...  
          x(node3)>=X_left_1 && x(node3)<=X_right_1 && ... 
          y(node1)>=Y_bottom_1 && y(node1)<=Y_top_1 && ...
          y(node2)>=Y_bottom_1 && y(node2)<=Y_top_1 && ...  
          y(node3)>=Y_bottom_1 && y(node3)<=Y_top_1    
               epsilon_r_final(e,1)=epsilon_r_1_dash-...
                   epsilon_r_1_dd*1i;
    elseif x(node1)>=X_left_2 && x(node1)<=X_right_2 && ...
          x(node2)>=X_left_2 && x(node2)<=X_right_2 && ...  
          x(node3)>=X_left_2 && x(node3)<=X_right_2 && ... 
          y(node1)>=Y_bottom_2 && y(node1)<=Y_top_2 && ...
          y(node2)>=Y_bottom_2 && y(node2)<=Y_top_2 && ...  
          y(node3)>=Y_bottom_2 && y(node3)<=Y_top_2   
                epsilon_r_final(e,1)=epsilon_r_2_dash-...
                   epsilon_r_2_dd*1i;
    elseif x(node1)>=X_left_3 && x(node1)<=X_right_3 && ...
           x(node2)>=X_left_3 && x(node2)<=X_right_3 && ...  
           x(node3)>=X_left_3 && x(node3)<=X_right_3 && ... 
           y(node1)>=Y_bottom_3 && y(node1)<=Y_top_3 && ...
           y(node2)>=Y_bottom_3 && y(node2)<=Y_top_3 && ...  
           y(node3)>=Y_bottom_3 && y(node3)<=Y_top_3    
                epsilon_r_final(e,1)=epsilon_r_3_dash-...
                   epsilon_r_3_dd*1i;
   elseif x(node1)>=X_left_4 && x(node1)<=X_right_4 && ...
           x(node2)>=X_left_4 && x(node2)<=X_right_4 && ...  
           x(node3)>=X_left_4 && x(node3)<=X_right_4 && ... 
           y(node1)>=Y_bottom_4 && y(node1)<=Y_top_4 && ...
           y(node2)>=Y_bottom_4 && y(node2)<=Y_top_4 && ...  
           y(node3)>=Y_bottom_4 && y(node3)<=Y_top_4    
                epsilon_r_final(e,1)=epsilon_r_4_dash-...
                   epsilon_r_4_dd*1i;
    elseif x(node1)>=X_left_5 && x(node1)<=X_right_5 && ...
           x(node2)>=X_left_5 && x(node2)<=X_right_5 && ...  
           x(node3)>=X_left_5 && x(node3)<=X_right_5 && ... 
           y(node1)>=Y_bottom_5 && y(node1)<=Y_top_5 && ...
           y(node2)>=Y_bottom_5 && y(node2)<=Y_top_5 && ...  
           y(node3)>=Y_bottom_5 && y(node3)<=Y_top_5    
               epsilon_r_final(e,1)=epsilon_r_5_dash-...
                   epsilon_r_5_dd*1i;
    else
        epsilon_r_final(e,1)=epsilon_r_air; 
    end
end 
tri_edge=[1 2;2 3;3 1];
local_edge= edge_array_local_2D(n,Ne,tri_edge);     
[global_edge,local_edge_array_elm_no ] = ...
        edge_array_global_2D_unique(local_edge);
No_of_edges=length(global_edge);  
edge_element = element_edge_array_2D(local_edge,global_edge,Ne); 
[ bc_edge_ABC_left bc_edge_ABC_right bc_ABC_left_edge_no bc_ABC_right_edge_no ]= ...
        find_port_boundary_edges(node_list,global_edge,X0,X1);
[ bc_edge_PEC_down bc_edge_PEC_top bc_PEC_down_edge_no bc_PEC_top_edge_no ]= ...
        find_PEC_boundary_edges(node_list,global_edge,Y0,Y1);    
% generating the edge sign array    
abs_edge=abs(edge_element);
sign_edge_element=zeros(Ne,3);
for ii=1:Ne
    for jj=1:3
        if edge_element(ii,jj)<0
            sign_edge_element(ii,jj)=-1;
        else 
            sign_edge_element(ii,jj)= 1;
         end
     end
end
% finding K_bar
[K_mean,f_mean]=obtain_K_mean_TM_mode ...
    (epsilon_r_final,freq,node_list,n,bc_node_ABC_left,bc_node_ABC_right, ...
    bc_node_ABC_down,bc_node_ABC_top,E0,edge_element,No_of_edges, ...
        bc_ABC_left_edge_no,bc_ABC_right_edge_no,bc_PEC_down_edge_no,bc_PEC_top_edge_no, ...
        sign_edge_element,tri_edge);
 
sigma_i_eps_r_dash=[0.1;0.1;0.1;0.1;0.1];
sigma_i_eps_r_dd=[0.1;0.1;0.1;0.1;0.1];  
% finding Ki
for ii=1:No_of_rv
    K_i(:,:,ii)=zeros(No_of_edges,No_of_edges);
end 
for e=1:Ne 
    b1=y(n(e,2))-y(n(e,3));
    b2=y(n(e,3))-y(n(e,1));
    b3=y(n(e,1))-y(n(e,2));
    b=[b1 b2 b3];
    c1=x(n(e,3))-x(n(e,2));
    c2=x(n(e,1))-x(n(e,3));
    c3=x(n(e,2))-x(n(e,1));
    c=[c1 c2 c3];
    Ae=0.5*(b1*c2-b2*c1);
    if Ae<0
        disp('negative area');
    end
    l1=sqrt((c3)^2+(b3)^2);
    l2=sqrt((c1)^2+(b1)^2);
    l3=sqrt((c2)^2+(b2)^2); 
    f11=(b1*b1)+(c1*c1);
    f12=(b1*b2)+(c1*c2);
    f13=(b1*b3)+(c1*c3);
    f21=(b2*b1)+(c2*c1);
    f22=(b2*b2)+(c2*c2);
    f23=(b2*b3)+(c3*c2);
    f31=(b3*b1)+(c3*c1);
    f32=(b3*b2)+(c2*c3);
    f33=(b3*b3)+(c3*c3);
    Fij=zeros(3,3);
    Fij(1,1)=(2*((l1)^2)*(f22-f12+f11))/(48*Ae);
    Fij(1,2)=(l1*l2*(f23-f22-(2*f13)+f12))/(48*Ae);
    Fij(1,3)=(l1*l3*(f21-(2*f23)-f11+f13))/(48*Ae);
    Fij(2,1)=(l1*l2*(f23- f22-(2*f13) +f12))/(48*Ae);
    Fij(2,2)=(2*((l2)^2)*(f33-f32+f22))/(48*Ae);
    Fij(2,3)=(l2*l3*(f31-f33-(2*f21)+f23))/(48*Ae);
    Fij(3,1)=(l1*l3*(f21-(2*f23)-f11+f13))/(48*Ae);
    Fij(3,2)=(l2*l3*(f31-f33-(2*f21)+f23))/(48*Ae);
    Fij(3,3)=(2*((l3)^2)*(f11-f13+f33))/(48*Ae);  
    node1=n(e,1);
    node2=n(e,2);
    node3=n(e,3); 
    for ii_1=1:No_of_rv
         Ke_temp(:,:,ii_1)=zeros(3,3);
    end
    if x(node1)>=X_left_1 && x(node1)<=X_right_1 && ...
          x(node2)>=X_left_1 && x(node2)<=X_right_1 && ...  
          x(node3)>=X_left_1 && x(node3)<=X_right_1 && ... 
          y(node1)>=Y_bottom_1 && y(node1)<=Y_top_1 && ...
          y(node2)>=Y_bottom_1 && y(node2)<=Y_top_1 && ...  
          y(node3)>=Y_bottom_1 && y(node3)<=Y_top_1    
           Ke_temp(:,:,1)= sigma_i_eps_r_dash(1,1)*Fij;
           Ke_temp(:,:,1+No_of_dielectric_layer)= -1i*sigma_i_eps_r_dd(1,1)*Fij;
    elseif x(node1)>=X_left_2 && x(node1)<=X_right_2 && ...
          x(node2)>=X_left_2 && x(node2)<=X_right_2 && ...  
          x(node3)>=X_left_2 && x(node3)<=X_right_2 && ... 
          y(node1)>=Y_bottom_2 && y(node1)<=Y_top_2 && ...
          y(node2)>=Y_bottom_2 && y(node2)<=Y_top_2 && ...  
          y(node3)>=Y_bottom_2 && y(node3)<=Y_top_2   
           Ke_temp(:,:,2)=sigma_i_eps_r_dash(2,1)*Fij;
           Ke_temp(:,:,2+No_of_dielectric_layer)= -1i*sigma_i_eps_r_dd(2,1)*Fij; 
    elseif x(node1)>=X_left_3 && x(node1)<=X_right_3 && ...
           x(node2)>=X_left_3 && x(node2)<=X_right_3 && ...  
           x(node3)>=X_left_3 && x(node3)<=X_right_3 && ... 
           y(node1)>=Y_bottom_3 && y(node1)<=Y_top_3 && ...
           y(node2)>=Y_bottom_3 && y(node2)<=Y_top_3 && ...  
           y(node3)>=Y_bottom_3 && y(node3)<=Y_top_3    
                Ke_temp(:,:,3)= sigma_i_eps_r_dash(3,1)*Fij;
           Ke_temp(:,:,3+No_of_dielectric_layer)= -1i*sigma_i_eps_r_dd(3,1)*Fij;
   elseif x(node1)>=X_left_4 && x(node1)<=X_right_4 && ...
           x(node2)>=X_left_4 && x(node2)<=X_right_4 && ...  
           x(node3)>=X_left_4 && x(node3)<=X_right_4 && ... 
           y(node1)>=Y_bottom_4 && y(node1)<=Y_top_4 && ...
           y(node2)>=Y_bottom_4 && y(node2)<=Y_top_4 && ...  
           y(node3)>=Y_bottom_4 && y(node3)<=Y_top_4    
                 Ke_temp(:,:,4)= sigma_i_eps_r_dash(4,1)*Fij;
           Ke_temp(:,:,4+No_of_dielectric_layer)= -1i*sigma_i_eps_r_dd(4,1)*Fij;
    elseif x(node1)>=X_left_5 && x(node1)<=X_right_5 && ...
           x(node2)>=X_left_5 && x(node2)<=X_right_5 && ...  
           x(node3)>=X_left_5 && x(node3)<=X_right_5 && ... 
           y(node1)>=Y_bottom_5 && y(node1)<=Y_top_5 && ...
           y(node2)>=Y_bottom_5 && y(node2)<=Y_top_5 && ...  
           y(node3)>=Y_bottom_5 && y(node3)<=Y_top_5    
              Ke_temp(:,:,5)=sigma_i_eps_r_dash(5,1)*Fij;
           Ke_temp(:,:,5+No_of_dielectric_layer)= -1i*sigma_i_eps_r_dd(5,1)*Fij;             
    end 
    for kk=1:No_of_rv
        K_matrix=K_i(:,:,kk);
        Ke_temp_matrix=-k0^2*Ke_temp(:,:,kk);
        for ii=1:3
          for jj=1:3 
                  K_matrix(abs_edge(e,ii),abs_edge(e,jj))=K_matrix(abs_edge(e,ii),abs_edge(e,jj))+ ...
                   (sign_edge_element(e,ii)*sign_edge_element(e,jj)* Ke_temp_matrix(ii,jj));
          end
        end
        K_i(:,:,kk)=K_matrix;
    end     
end 
% Applying Boundary condition for K1 and K2  
% Applying PEC boundary conditions
bc_PEC_edge_no_concate=[ bc_PEC_down_edge_no ; bc_PEC_top_edge_no];
bc_PEC_edge_no=unique(bc_PEC_edge_no_concate);
for ii=1:No_of_rv
    K_matrix=K_i(:,:,ii);
    K_matrix(bc_PEC_edge_no ,:)=[]; % Deleting corresponding row 
    K_matrix(:,bc_PEC_edge_no)=[]; % Deleting corresponding column
    K_i_bc(:,:,ii)=K_matrix;
end
monte_carlo_iteration=10000;
gaussian_RV=zeros(monte_carlo_iteration,No_of_dielectric_layer);
for ii=1:No_of_rv
    seed_value=ii;
    gaussian_RV(:,ii)=random_number_seed_gaussian(monte_carlo_iteration,seed_value);
end
xi_value = gaussian_RV(:,1:No_of_rv); 
 No_of_LHS_samples=25;
%No_of_LHS_samples=10;
[lhs_samples_2,min_length_lhs_2]= ...
        obtain_lhs_samples_final(No_of_LHS_samples ,No_of_rv ,xi_value);
% obtaining stochastic response for k samples 
U_matrix=zeros(No_of_edges-length(bc_PEC_edge_no),min_length_lhs_2);
for kk=1:min_length_lhs_2
   kk
   epsilon_r_final=zeros(Ne,1);
   lhs_samples_itr=lhs_samples_2(kk,:);
   for e=1:Ne 
       node1=n(e,1);
        node2=n(e,2);
        node3=n(e,3); 
        if x(node1)>=X_left_1 && x(node1)<=X_right_1 && ...
              x(node2)>=X_left_1 && x(node2)<=X_right_1 && ...  
              x(node3)>=X_left_1 && x(node3)<=X_right_1 && ... 
              y(node1)>=Y_bottom_1 && y(node1)<=Y_top_1 && ...
              y(node2)>=Y_bottom_1 && y(node2)<=Y_top_1 && ...  
              y(node3)>=Y_bottom_1 && y(node3)<=Y_top_1    
                   epsilon_r_final(e,1)=epsilon_r_1_dash ...
                          +sigma_i_eps_r_dash(1,1)*lhs_samples_itr(1,1)-...
                          (1i*(epsilon_r_1_dd+sigma_i_eps_r_dd(1,1)* ...
                          lhs_samples_itr(1,1+No_of_dielectric_layer))); 
        elseif x(node1)>=X_left_2 && x(node1)<=X_right_2 && ...
              x(node2)>=X_left_2 && x(node2)<=X_right_2 && ...  
              x(node3)>=X_left_2 && x(node3)<=X_right_2 && ... 
              y(node1)>=Y_bottom_2 && y(node1)<=Y_top_2 && ...
              y(node2)>=Y_bottom_2 && y(node2)<=Y_top_2 && ...  
              y(node3)>=Y_bottom_2 && y(node3)<=Y_top_2   
                    epsilon_r_final(e,1)=epsilon_r_2_dash ...
                          +sigma_i_eps_r_dash(2,1)*lhs_samples_itr(1,2)-...
                          (1i*(epsilon_r_2_dd+sigma_i_eps_r_dd(2,1)* ...
                          lhs_samples_itr(1,2+No_of_dielectric_layer))); 
        elseif x(node1)>=X_left_3 && x(node1)<=X_right_3 && ...
               x(node2)>=X_left_3 && x(node2)<=X_right_3 && ...  
               x(node3)>=X_left_3 && x(node3)<=X_right_3 && ... 
               y(node1)>=Y_bottom_3 && y(node1)<=Y_top_3 && ...
               y(node2)>=Y_bottom_3 && y(node2)<=Y_top_3 && ...  
               y(node3)>=Y_bottom_3 && y(node3)<=Y_top_3    
                    epsilon_r_final(e,1)=epsilon_r_3_dash ...
                          +sigma_i_eps_r_dash(3,1)*lhs_samples_itr(1,3)-...
                          (1i*(epsilon_r_3_dd+sigma_i_eps_r_dd(3,1)* ...
                          lhs_samples_itr(1,3+No_of_dielectric_layer)));  
       elseif x(node1)>=X_left_4 && x(node1)<=X_right_4 && ...
               x(node2)>=X_left_4 && x(node2)<=X_right_4 && ...  
               x(node3)>=X_left_4 && x(node3)<=X_right_4 && ... 
               y(node1)>=Y_bottom_4 && y(node1)<=Y_top_4 && ...
               y(node2)>=Y_bottom_4 && y(node2)<=Y_top_4 && ...  
               y(node3)>=Y_bottom_4 && y(node3)<=Y_top_4    
                    epsilon_r_final(e,1)=epsilon_r_4_dash ...
                          +sigma_i_eps_r_dash(4,1)*lhs_samples_itr(1,4)-...
                          (1i*(epsilon_r_4_dd+sigma_i_eps_r_dd(4,1)* ...
                          lhs_samples_itr(1,4+No_of_dielectric_layer)));  
        elseif x(node1)>=X_left_5 && x(node1)<=X_right_5 && ...
               x(node2)>=X_left_5 && x(node2)<=X_right_5 && ...  
               x(node3)>=X_left_5 && x(node3)<=X_right_5 && ... 
               y(node1)>=Y_bottom_5 && y(node1)<=Y_top_5 && ...
               y(node2)>=Y_bottom_5 && y(node2)<=Y_top_5 && ...  
               y(node3)>=Y_bottom_5 && y(node3)<=Y_top_5    
                   epsilon_r_final(e,1)=epsilon_r_5_dash ...
                          +sigma_i_eps_r_dash(5,1)*lhs_samples_itr(1,5)-...
                          (1i*(epsilon_r_5_dd+sigma_i_eps_r_dd(5,1)* ...
                          lhs_samples_itr(1,5+No_of_dielectric_layer)));  
        else
            epsilon_r_final(e,1)=epsilon_r_air; 
        end       
    end 
    [Kij,f_mean]=obtain_K_mean_TM_mode ...
        (epsilon_r_final,freq,node_list,n,bc_node_ABC_left,bc_node_ABC_right, ...
        bc_node_ABC_down,bc_node_ABC_top,E0,edge_element,No_of_edges, ...
            bc_ABC_left_edge_no,bc_ABC_right_edge_no,bc_PEC_down_edge_no,bc_PEC_top_edge_no, ...
            sign_edge_element,tri_edge);
    E_field_samples =Kij\f_mean;         
    U_matrix(:,kk)=E_field_samples;
end
[SVD_U,SVD_Sigma,SVD_V_T]=svd(U_matrix);
sing_values=diag(SVD_Sigma);
count_sing_value=1;
N_value=length(sing_values);

% for ii=1:N_value
%      RIC(ii,1)=sum(sing_values(ii+1:N_value).*sing_values(ii+1:N_value))/sum(sing_values.*sing_values);
%      if RIC(ii,1) <=1/100000 && break_loop==0
%          m_count=ii;
%          break_loop=1;
%      end
% end
for ii=1:N_value
     RIC_value(ii,1)=sum(sing_values(1:ii).*sing_values(1:ii))/sum(sing_values.*sing_values);
end

err_value_1=0.9999999;
err_value_2=0.9999990;
err_value_3=0.9999900;
err_value_4=0.9999000;
err_value_5=0.9990000;
err_value_6=0.9900000;
err_value_7=0.9000000;
err_value_8=0.8000000;
err_value=[err_value_1;err_value_2;err_value_3;err_value_4;...
    err_value_5;err_value_6;err_value_7;err_value_8];
break_loop=0;
for ii=1:N_value 
     if RIC_value(ii,1) >=err_value_1 && break_loop==0
         m_count=ii;
         break_loop=1;
     end
end
m_count_list=zeros(length(err_value),1);
for jj=1:length(err_value)
    break_loop=0;
    for ii=1:N_value 
         if RIC_value(ii,1) >=err_value(jj,1) && break_loop==0
             m_count_list(jj,1)=ii;
             break_loop=1;
         end
    end
end
U_POD=SVD_U(:,1:m_count);
K_bar=transpose(U_POD)*K_mean*U_POD;
f_bar=transpose(U_POD)*f_mean;
for kk=1:No_of_rv
    K_i_bar(:,:,kk)=transpose(U_POD)*K_i_bc(:,:,kk)*U_POD;
end
lower_limit=Y0/100;
upper_limit=Y1/100;
left_Ey_perturb=zeros(monte_carlo_iteration,1);
right_Ey_perturb=zeros(monte_carlo_iteration,1);
for kk=1:monte_carlo_iteration  
    kk
    K_final_bar=K_bar;
    for jj=1:No_of_rv
        K_final_bar=K_final_bar+K_i_bar(:,:,jj)*xi_value(kk,jj);
    end 
    E_field_POD_coefficient=K_final_bar\f_bar;
    E_field_final_POD=U_POD*E_field_POD_coefficient;
    % Adding the unknown quantity obtained with the known quantity
    count=1;
    edge_count=1;
    E_field_final=zeros(No_of_edges,1);
     for ii=1:No_of_edges
        if edge_count<=length(bc_PEC_edge_no)
            if(ii==bc_PEC_edge_no(edge_count) )
                E_field_final(ii)=0; 
                edge_count=edge_count+1;
            else
                E_field_final(ii)=E_field_final_POD(count);
                count=count+1; 
            end 
        else
            E_field_final(ii)=E_field_final_POD(count);
            count=count+1;        
        end
     end
       
    Ey_final_port1=E_field_final(bc_ABC_left_edge_no);
    Ey_final_port2=E_field_final(bc_ABC_right_edge_no);
    [sign_value ]=get_sign_element_port(global_edge(bc_ABC_left_edge_no,:),n(element_left_bc,:));
    Ey_final_port1=Ey_final_port1.*sign_value;
    [sign_value ]=get_sign_element_port(global_edge(bc_ABC_right_edge_no,:),n(element_right_bc,:));
    Ey_final_port2=Ey_final_port2.*sign_value;
    y=node_list(:,2)/100;
    y_node_value_left=(y(global_edge(bc_ABC_left_edge_no,1))+ ...
            y(global_edge(bc_ABC_left_edge_no,2)))/2;
    N_value_left=length(bc_ABC_left_edge_no);
    y_node_value_right=(y(global_edge(bc_ABC_right_edge_no,1))+ ...
            y(global_edge(bc_ABC_right_edge_no,2)))/2;
    N_value_right=length(bc_ABC_right_edge_no);
    Ey_final_port1_function= ...
        interpolation_lagrange(Ey_final_port1,N_value_left,y_node_value_left);
    Ey_final_port2_function= ...
        interpolation_lagrange(Ey_final_port2,N_value_right,y_node_value_right);
    function_k_1 = @(y) Ey_final_port1_function(y)*E0; 
    left_Ey_perturb(kk,1)=perform_1D_integration(function_k_1,lower_limit,upper_limit); 
    function_p_1=@(y) Ey_final_port2_function(y)*E0*exp(-i*k0*X1/100);  
    right_Ey_perturb(kk,1)=perform_1D_integration(function_p_1,lower_limit,upper_limit);   
end
function_k_2=@(y) (E0^2);
gamma_fem=(left_Ey_perturb-perform_1D_integration(function_k_2,lower_limit,upper_limit)) ...
          /perform_1D_integration(function_k_2,lower_limit,upper_limit);
mag_gamma_fem_POD=abs(gamma_fem);
function_p_2=@(y) (E0^2)*exp(-i*k0*2*X1/100);
trans_fem=right_Ey_perturb/perform_1D_integration(function_p_2,lower_limit,upper_limit);
trans_fem_abs_POD =abs(trans_fem);  
time_elapsed_POD_2D_multilayer=toc;
[prob_dist_vector_fem_POD, set_of_points_fem_POD ]=ksdensity(mag_gamma_fem_POD);
figure(1)
colorstring = 'krbgmwy'; 
plot(set_of_points_fem_POD,prob_dist_vector_fem_POD ,'LineWidth',1.5,'Color',colorstring(3));
xlabel('Magnitude of Reflection Coefficient  , $|\Gamma| $ ','fontsize',16,'Interpreter','latex');
ylabel('Probability density function  ','fontsize',16,'Interpreter','latex'); 

[prob_dist_vector_fem_trans_coeff_POD, set_of_points_fem_trans_coeff_POD ]=ksdensity(trans_fem_abs_POD);
figure(2)
colorstring ='krbgmwy';
plot(set_of_points_fem_trans_coeff_POD,prob_dist_vector_fem_trans_coeff_POD , 'LineWidth',1.5,'Color',colorstring(3));
xlabel('Magnitude of Tranmission Coefficient  , $|T| $ ','fontsize',16,'Interpreter','latex');
ylabel('Probability density function  ','fontsize',16,'Interpreter','latex'); 
