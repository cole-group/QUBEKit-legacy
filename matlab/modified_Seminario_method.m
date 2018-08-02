function modified_Seminario_method( inputfilefolder,  outputfilefolder, vibrational_scaling)

%   Program to implement the Modified Seminario Method
%   Written by Alice E. A. Allen, TCM, University of Cambridge
%   Reference using ... 

% Pass as arguements the input folder containing the zmat.log/lig.log,
% lig.fchk and optional Zmat.z file, the output folder where the new
% parameters will be written and the vibrational frequency scaling constant
% required. 

%Create log file
fid_log = fopen(horzcat(outputfilefolder,'MSM_log'), 'wt');
fprintf(fid_log, '%s\n', 'Modified Seminario Method ');
fprintf(fid_log, '%s %s\n', 'Parametrization started for files in folder', inputfilefolder);
fprintf(fid_log, '%s %s\n', 'Time is now: ', datestr(clock, 0));

%Square the vibrational scaling used for frequencies
vibrational_scaling_squared = vibrational_scaling^2; 

%Import all input data
[ bond_list,angle_list, coords,N, hessian,atom_names ] = input_data_processing( inputfilefolder );

%Find bond lengths
bond_lengths  = zeros(N,N);

for i = 1:N
    for j = 1:N
        bond_lengths(i,j) = norm(coords(:,i) - coords(:,j));
    end
end

%Eigenvectors and eigenvalues calculated
eigenvectors = zeros(3,3,N,N);
eigenvalues = zeros(N,N,3);

for i = 1:N
    for j = 1:N
        [A,B] = eig(hessian( ((i - 1)*3 + 1):((i)*3), ( (j - 1)*3 + 1):((j)*3) ) );
        eigenvalues(i,j,:) = diag(B);
        eigenvectors(:,:,i,j) = (A);
    end
end

% The bond values are calculated and written to file
[unique_values_bonds] = bonds_calculated_printed( outputfilefolder, vibrational_scaling_squared, bond_list, bond_lengths, atom_names, eigenvalues, eigenvectors, coords );

% The angle values are calculated and written to file
[unique_values_angles] = angles_calculated_printed( outputfilefolder, vibrational_scaling_squared, angle_list, bond_lengths, atom_names, eigenvalues, eigenvectors, coords );

%The final section finds the average bond and angle terms for each
%bond/angle class if the .z exists to supply angle/bond classes and then 
%writes the new terms to a .sb file
if exist( horzcat(inputfilefolder,'Zmat.z'), 'file') == 2
   
    average_values_across_classes( unique_values_bonds, unique_values_angles, outputfilefolder )
   
    %The functions below creates a new .sb file with the modified Seminario
    %terms
    sb_file_new_parameters(outputfilefolder, 'Modified_Scaled');
    
end

%Log file write
fprintf(fid_log, '%s %s\n', 'Parametrization finished for files in folder', inputfilefolder);
fprintf(fid_log, '%s %s\n', 'Time is now: ', datestr(clock, 0));
fclose(fid_log);

end

