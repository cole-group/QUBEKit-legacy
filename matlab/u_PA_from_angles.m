function [ u_PA ] = u_PA_from_angles( atom_A, atom_B, atom_C, coords )
%This gives the vector in the plane A,B,C and perpendicular to A to B

diff_AB = coords(:,atom_B) - coords(:, atom_A);
u_AB = (diff_AB) ./ norm(diff_AB); 

diff_CB = coords(:,atom_B) - coords(:, atom_C);
u_CB = (diff_CB) ./ norm(diff_CB); 

u_N = unit_vector_N( u_CB, u_AB );

u_PA = cross(u_N,  u_AB);
u_PA = u_PA /norm(u_PA);

end

