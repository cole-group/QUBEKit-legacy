function [ unique_values_angles ] = angles_calculated_printed( outputfilefolder, vibrational_scaling_squared, angle_list, bond_lengths, atom_names, eigenvalues, eigenvectors, coords )
%This function uses the modified Seminario method to find the angle
%parameters and print them to file

%Open output file angle parameters are written to 
fid = fopen(horzcat(outputfilefolder,'Modified_Seminario_Angle'), 'w');

k_theta = zeros(1, size(angle_list,1));
theta_0 = zeros(1, size(angle_list,1));
unique_values_angles = cell(6,size(angle_list,1));


%Connectivity information for Modified Seminario Method
central_atoms_angles = cell(max(max(angle_list)));

%A structure is created with the index giving the central atom of the
%angle, an array then lists the angles with that central atom. 
%ie. central_atoms_angles{3} contains an array of angles with central atom
%3
for i =1:size(angle_list,1)
    %For angle ABC, atoms A C are written to array
    central_atoms_angles{angle_list(i,2)}( end  + 1, 1) = angle_list(i, 1);
    central_atoms_angles{angle_list(i,2)}(end, 2) =  angle_list(i, 3);
    central_atoms_angles{angle_list(i,2)}(end, 3) =  i; % Position in angle list tracked 
    
    %For angle ABC, atoms C A are written to array
    central_atoms_angles{angle_list(i,2)}( end  + 1, 1) = angle_list(i, 3);
    central_atoms_angles{angle_list(i,2)}(end, 2) =  angle_list(i, 1);
    central_atoms_angles{angle_list(i,2)}(end, 3) =  i;
end

%Sort rows by atom numbers
for i = 1:size(central_atoms_angles,1)
    central_atoms_angles{i} = sortrows(central_atoms_angles{i});
end

%Find normals u_PA for each angle
unit_PA_all_angles = cell(max(max(angle_list)));
for i = 1:size(central_atoms_angles,1)
    for j =1:size( central_atoms_angles{i},1 )
        %For the angle at central_atoms_angles{i}(j,:) the corresponding
        %u_PA value is found for the plane ABC and bond AB, where ABC
        %corresponds to the order of the arguements
        %This is why the reverse order was also added
        unit_PA_all_angles{i}(j,:) = u_PA_from_angles(central_atoms_angles{i}(j,1), i, central_atoms_angles{i}(j,2), coords);
    end
end

%Finds the contributing factors from the other angle terms
scaling_factor_all_angles = cell(max(max(angle_list))); %This array will contain scaling factor and angle list position
for i = 1:size(central_atoms_angles,1)
    for j =1:size( central_atoms_angles{i},1 )
        n = 1; 
        m = 1; 
        angles_around = 0; 
        additional_contributions = 0; 
        
        %Position in angle list
        scaling_factor_all_angles{i}(j,2) =  central_atoms_angles{i}(j,3); 
        
        % Goes through the list of angles with the same central atom
        % And computes the term need for the modified Seminario method
        
        %Forwards directions, finds the same bonds with the central atom i  
        while( ( (j + n ) <= size(central_atoms_angles{i},1) )&&central_atoms_angles{i}(j,1) == central_atoms_angles{i}(j+n,1) )
            additional_contributions = additional_contributions + (abs(dot(unit_PA_all_angles{i}(j,:), unit_PA_all_angles{i}((j + n),:)))).^2 ;            
            n = n + 1; 
            angles_around =  angles_around + 1; 
        end
        
        %Backwards direction, finds the same bonds with the central atom i   
        while( ( (j - m ) >= 1 )&& central_atoms_angles{i}(j,1) == central_atoms_angles{i}(j-m,1) )
            additional_contributions = additional_contributions + (abs(dot(unit_PA_all_angles{i}(j,:), unit_PA_all_angles{i}((j - m),:)))).^2;            
            m = m + 1; 
            angles_around =  angles_around + 1; 
        end
        
        if (n ~= 1 || m ~= 1)
            %Finds the mean value of the additional contribution
            %To change to normal Seminario method comment out + part 
            scaling_factor_all_angles{i}(j,1) = 1 + (additional_contributions)/((m  + n - 2) );  
        else 
            scaling_factor_all_angles{i}(j,1) = 1;
        end
    end
end

%Orders the scaling factors according to the angle list
scaling_factors_angles_list =  cell(size(angle_list,1));  
for i = 1:size(central_atoms_angles,1)
    for j =1:size( central_atoms_angles{i},1 )
        scaling_factors_angles_list{scaling_factor_all_angles{i}(j,2)}(end + 1) = scaling_factor_all_angles{i}(j,1); 
    end
end



%Finds the angle force constants with the scaling factors included for each
%angle
for i = 1:size(angle_list,1)
    %Ensures that there is no difference when the ordering is changed 
    [AB_k_theta, AB_theta_0] = force_angle_constant( angle_list(i,1), angle_list(i,2), angle_list(i,3), bond_lengths, eigenvalues, eigenvectors, coords, scaling_factors_angles_list{i}(1), scaling_factors_angles_list{i}(2) ); 
    [BA_k_theta, BA_theta_0] = force_angle_constant( angle_list(i,3), angle_list(i,2), angle_list(i,1), bond_lengths, eigenvalues, eigenvectors, coords, scaling_factors_angles_list{i}(2), scaling_factors_angles_list{i}(1) ); 
    k_theta(i) = ( AB_k_theta + BA_k_theta ) / 2;
    
    %Vibrational_scaling takes into account DFT deficities/ anharmocity 
    k_theta(i) =  k_theta(i) * vibrational_scaling_squared;
    
    theta_0(i) = ( AB_theta_0 +  BA_theta_0 ) / 2;
        
    fprintf(fid, '%s%s%s%s%s    % 12.3f      % 12.3f     %i %i %i\n', char(atom_names{angle_list(i,1)}), '-', char(atom_names{angle_list(i,2) }),'-', char(atom_names{angle_list(i,3) }), k_theta(i), theta_0(i), angle_list(i,1), angle_list(i,2), angle_list(i,3));
    
    unique_values_angles{1,i} = char(atom_names{angle_list(i,1)});
    unique_values_angles{2,i} = char(atom_names{angle_list(i,2)});
    unique_values_angles{3,i} = char(atom_names{angle_list(i,3)});
    unique_values_angles{4,i} = k_theta(i);
    unique_values_angles{5,i} = theta_0(i);
    unique_values_angles{6,i} = 1;
end

end

