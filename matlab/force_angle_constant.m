function [ k_theta, theta_0 ] = force_angle_constant( atom_A, atom_B, atom_C, bond_lengths, eigenvalues, eigenvectors, coords, scaling_1, scaling_2 )
%Force Constant- Equation 14 of seminario calculation paper - gives force
%constant for angle (in kcal/mol/rad^2) and equilibrium angle in degrees

%Vectors along bonds calculated
diff_AB = coords(:,atom_B) - coords(:, atom_A);
u_AB = (diff_AB) ./ norm(diff_AB); 

diff_CB = coords(:,atom_B) - coords(:, atom_C);
u_CB = (diff_CB) ./ norm(diff_CB); 

%Bond lengths and eigenvalues found
bond_length_AB = bond_lengths(atom_A,atom_B);
eigenvalues_AB = eigenvalues(atom_A, atom_B, :);
eigenvectors_AB = eigenvectors(1:3, 1:3, atom_A, atom_B);

bond_length_BC = bond_lengths(atom_B,atom_C);
eigenvalues_CB = eigenvalues(atom_C, atom_B,  :);
eigenvectors_CB = eigenvectors(1:3, 1:3, atom_C, atom_B);

%Normal vector to angle plane found
u_N = unit_vector_N( u_CB, u_AB );

u_PA = cross(u_N,  u_AB);
u_PA = u_PA /norm(u_PA);

u_PC = cross(u_CB, u_N);
u_PC = u_PC /norm(u_PC);

sum_first = zeros(1,3); 
sum_second = zeros(1,3); 


%Projections of eigenvalues
for i = 1:3
    eig_AB_i = (eigenvectors_AB(:,i));
    sum_first(i) = ((eigenvalues_AB(i)) * abs( (dot(u_PA, eig_AB_i))) ); 
    eig_BC_i = (eigenvectors_CB(:,i));
    sum_second(i) =  ((eigenvalues_CB(i)) * abs(dot(u_PC, eig_BC_i))) ; 
end

%Sum contributions
sum_first = sum((sum_first)); 
sum_second = sum((sum_second)); 
 
%Scaling due to additional angles - Modified Seminario Part
sum_first = sum_first/scaling_1; 
sum_second = sum_second/scaling_2; 

%Added as two springs in series
k_theta = (( 1 / ( bond_length_AB^2 * sum_first) ) + ( 1 / (bond_length_BC^2 * sum_second) ) );  
k_theta = 1/(k_theta);

k_theta = - k_theta; %Change to OPLS form

k_theta = abs(k_theta  * 0.5); %Change to OPLS form

%Equilibrium Angle
theta_0 = acosd(dot(u_AB, u_CB));

%If the vectors u_CB and u_AB are linearly dependent u_N cannot be defined
% This case is dealt with here
if isnan( sum(u_N) )
    scaling_1 = 1; 
    scaling_2 = 1; 
    [ k_theta, theta_0 ] = force_angle_constant_special_case( atom_A, atom_B, atom_C, bond_lengths, eigenvalues, eigenvectors, coords, scaling_1, scaling_2 );
end

end

