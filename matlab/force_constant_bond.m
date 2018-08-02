function [ k_AB ] = force_constant_bond(atom_A, atom_B,eigenvalues, eigenvectors, coords )
%Force Constant - Equation 10 of Seminario paper - gives force
%constant for bond

%Eigenvalues and eigenvectors calculated
eigenvalues_AB = eigenvalues(atom_A,atom_B,  :);
eigenvectors_AB = eigenvectors(1:3, 1:3, atom_A, atom_B  );

%Vector along bond
diff_AB = coords(:,atom_B) - coords(:, atom_A);
unit_vectors_AB = (diff_AB) / norm(diff_AB); 

k_AB = 0; 

%Projection of eigenvalues
for i = 1:3
    k_AB = k_AB + (eigenvalues_AB(i)) .* abs(dot(unit_vectors_AB, eigenvectors_AB(:,i))) ; 
end

k_AB = - k_AB ; %Convert to OPLS form  
k_AB = k_AB * 0.5; %Convert to OPLS form

end

