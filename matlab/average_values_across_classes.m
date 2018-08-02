function average_values_across_classes( unique_values_bonds, unique_values_angles, outputfilefolder )
%This function finds the average bond and angle term for each class. 

 remove_rows = [];
    
    %Find Average Values Bonds
    for i = 1:(size(unique_values_bonds,2) - 1)
        for j = (i+1):size(unique_values_bonds,2)
            %Finds if the bond class has already been encountered
            if ( strcmp(unique_values_bonds{1,i}, unique_values_bonds(1,j)) && strcmp(unique_values_bonds{2,i}, unique_values_bonds(2,j) ) ) || ( strcmp(unique_values_bonds{1,i}, unique_values_bonds(2,j)) && strcmp(unique_values_bonds{2,i}, unique_values_bonds(1,j) ) )
                unique_values_bonds{3,i} = unique_values_bonds{3,i} + unique_values_bonds{3,j};
                unique_values_bonds{4,i} = unique_values_bonds{4,i} + unique_values_bonds{4,j};
                unique_values_bonds{5,i} = unique_values_bonds{5,i} + 1;
                remove_rows(end+1) = j;
            end
        end
    end
    
    %Removes bonds classes that were already present
    unique_values_bonds(:,unique(remove_rows)) = [];
    
    % Finds mean value 
    for i= 1:size(unique_values_bonds,2)
        unique_values_bonds{3,i} = unique_values_bonds{3,i} / unique_values_bonds{5,i};
        unique_values_bonds{4,i} = unique_values_bonds{4,i} / unique_values_bonds{5,i};
    end
    
    %Average Bonds Printed  
    fid = fopen(horzcat(outputfilefolder,'Average_Modified_Seminario_Bonds'), 'w');
    
    for i = 1:size(unique_values_bonds,2)
        fprintf(fid, '%s%s%s,  %-8.2f,      %-8.3f  \n', unique_values_bonds{1,i}, '-', unique_values_bonds{2,i}, (unique_values_bonds{3,i}), (unique_values_bonds{4,i}) );
    end
    
    fclose(fid);
    
    %Find Average Values Angles
    remove_rows_angles = [];
    
    for i = 1:(size(unique_values_angles,2) - 1)
        for j = (i+1):size(unique_values_angles,2)
            %Finds if the angle class has already been encountered
            if ( strcmp(unique_values_angles{1,i}, unique_values_angles(1,j)) && strcmp(unique_values_angles{2,i}, unique_values_angles(2,j) )  && strcmp(unique_values_angles{3,i}, unique_values_angles(3,j) ) ) || ( strcmp(unique_values_angles{1,i}, unique_values_angles(3,j)) && strcmp(unique_values_angles{2,i}, unique_values_angles(2,j) ) && strcmp(unique_values_angles{3,i}, unique_values_angles(1,j) ) )
                unique_values_angles{4,i} = unique_values_angles{4,i} + unique_values_angles{4,j};
                unique_values_angles{5,i} = unique_values_angles{5,i} + unique_values_angles{5,j};
                unique_values_angles{6,i} = unique_values_angles{6,i} + 1;
                remove_rows_angles(end+1) = j;
            end
        end
    end
    
    %Removes angle classes that were already present
    unique_values_angles(:,unique(remove_rows_angles)) = [];
    
    %Find the mean angle values
    for i= 1:size(unique_values_angles,2)
        unique_values_angles{4,i} = unique_values_angles{4,i} / unique_values_angles{6,i};
        unique_values_angles{5,i} = unique_values_angles{5,i} / unique_values_angles{6,i};
    end
    
    %Average Angles Printed   
    fid = fopen(horzcat(outputfilefolder,'Average_Modified_Seminario_Angles'), 'wt');
    
    for i = 1:size(unique_values_angles,2)
        fprintf(fid, '%s%s%s%s%s, %-12.3f,      %-12.3f  \n', unique_values_angles{1,i}, '-', unique_values_angles{2,i}, '-', unique_values_angles{3,i}, unique_values_angles{4,i} , unique_values_angles{5,i} );
    end
    
    fclose(fid);
    
end

