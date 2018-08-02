function [ unique_values_bonds ] = bonds_calculated_printed( outputfilefolder, vibrational_scaling_squared, bond_list, bond_lengths, atom_names, eigenvalues, eigenvectors, coords )
%This function uses the Seminario method to find the bond
%parameters and print them to file

%Open output file bond parameters are written to 
fid = fopen(horzcat(outputfilefolder,'Modified_Seminario_Bonds'), 'w');

k_b = zeros(1, size(bond_list,1));
bond_length_list = zeros(1, size(bond_list,1));
unique_values_bonds = cell(5,size(bond_list,1)); %used to find average values 

for i = 1:size(bond_list,1)
    AB = force_constant_bond(bond_list(i,1), bond_list(i,2),eigenvalues, eigenvectors, coords );
    BA = force_constant_bond(bond_list(i,2), bond_list(i,1),eigenvalues, eigenvectors, coords );
    k_b(i) = ( AB + BA ) /2; % Order of bonds sometimes causes slight differences, find the mean
    
    %Vibrational_scaling takes into account DFT deficities/ anharmocity    
    k_b(i) =  k_b(i) * vibrational_scaling_squared;
    
    bond_length_list(i) =  bond_lengths(bond_list(i,1), bond_list(i,2) );
    fprintf(fid, '%s%s%s %8.3f      %8.3f      %i %i \n', char(atom_names{bond_list(i,1)}), '-', char(atom_names{bond_list(i,2) }), k_b(i), (bond_length_list(i)), bond_list(i,1), bond_list(i,2));
    
    unique_values_bonds{1,i} = char(atom_names{bond_list(i,1)});
    unique_values_bonds{2,i} = char(atom_names{bond_list(i,2)});
    unique_values_bonds{3,i} = k_b(i);
    unique_values_bonds{4,i} = bond_length_list(i);
    unique_values_bonds{5,i} = 1;
end

fclose(fid);

end

