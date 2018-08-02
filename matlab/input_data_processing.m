function [ bond_list,angle_list, coords,N, hessian,atom_names ] = input_data_processing( inputfilefolder )
%This function takes input data that is need from the files supplied
%Function extracts input coords and hessian from .fchk file, bond and angle
%lists from .log file and atom names if a z-matrix is supplied

%Gets Hessian in unprocessed form and writes .xyz file too 
[unprocessed_Hessian,N] =  coords_from_fchk( inputfilefolder, 'lig.fchk' );

delimiterIn = ' ';
headerlinesIn = 2;
coords = importdata(horzcat(inputfilefolder,'input_coords.xyz'), delimiterIn, headerlinesIn);

%Get bond and angle lists
[bond_list,angle_list] =  bond_angle_lists(inputfilefolder);

%Get  OPLS number to name list and format correctly
OPLS_number_to_name = importdata('./Number_to_Atom_type');
OPLS_number_to_name = strtrim(OPLS_number_to_name);

for i = 1:size(OPLS_number_to_name,1)
    OPLS_number_to_name{i} = strsplit(OPLS_number_to_name{i});
end

%Imports Zmat and extract just Z matrix part
if exist( horzcat(inputfilefolder,'Zmat.z'), 'file') == 2
    fid = fopen(horzcat(inputfilefolder,'Zmat.z')); %Boss types Zmat is fine
    fidalone = fopen(horzcat(inputfilefolder,'Zmat_alone'), 'wt'); % Script makes this file
    fid_atom_names = fopen(horzcat(inputfilefolder,'Zmat_names'), 'wt'); % Script makes this file
    
    
    tline = fgets(fid);
    tline = fgets(fid);
    
    %Find number of dummy atoms
    number_dummy = 0;
    tmp = strsplit(strtrim(tline));
    
    while strcmp(tmp{3},  '-1')
        number_dummy = number_dummy +1;
        tline = fgets(fid);
        tmp = strsplit(strtrim(tline));
    end
    
    for i=1:(N)
        fprintf(fidalone, '%s', tline);
        tline = fgets(fid);
    end
    
    fclose(fidalone);
    
    print_out = 'n';
    while ischar(tline)
        tline = fgets(fid);
        
        if size(tline,2) > 5 && strcmp(tline(1:6), ' Final')
            tline = fgets(fid);
            tline = fgets(fid);
            print_out = 'y';
        end
        
        if strcmp(print_out, 'y')
            if ischar(tline)
                fprintf(fid_atom_names, '%s', tline);
            end
        end
    end
    
    fclose(fid);
    fclose(fid_atom_names);
    
    zmat= importdata(horzcat(inputfilefolder,'Zmat_alone')); 
end

%Process Hessian to full matrix
length_hessian = 3 * N;

hessian = zeros(length_hessian, length_hessian); %In proper format

unprocessed_Hessian = (unprocessed_Hessian .* (627.509391) )./ (0.529.^2) ; %Change from Hartree/bohr to kcal/mol /ang

for i = 1:length_hessian
    hessian(i,1:i) = unprocessed_Hessian(( 0.5 * i * (i - 1) + 1) : ( 0.5 * (i )* ( i + 1) ) );
    hessian(1:i, i) = unprocessed_Hessian(( 0.5 * i * (i - 1) + 1) : ( 0.5 * (i )* ( i + 1) ) );
end

%xyz Input coordinates processed
if size(coords.data,2) == 3
    atom_element = coords.textdata;
    atom_element(1:2,: )  = [];
    coords = coords.data;
    coords = coords';
else
    if size(coords.data,2) == 4
        atom_element = coords.textdata;
        coords = coords.data;
        coords = coords';
        coords  = coords(2:4, :);
    else
        fprintf('Problem with xyz file');
    end
end


%Atoms names from Zmat
if exist( horzcat(inputfilefolder,'Zmat.z'), 'file') == 2
    
    if size(zmat.textdata,2) > 2
        atom_numbers = zmat.textdata(:,3);
    else
        atom_numbers = zmat.data(:,1);
    end
    
    %Ensure atom_numbers are the correct type
    if  iscell(atom_numbers)
        atom_numbers = str2num(char(atom_numbers));
    end
    
    atom_names = zmat.textdata(:,2);
    
    if max(atom_numbers) < 800 %Get new names
        for i = 1:size(atom_names,1)
            for j = 1:size(OPLS_number_to_name,1)
                if str2double(char(OPLS_number_to_name{j}(1))) == (atom_numbers(i))
                    atom_names{i}  = char(OPLS_number_to_name{j}(2));
                end
            end
        end
    else
        %If CM1A style Z matrix then takes the names from the end of the
        %file
        zmat_names= importdata(horzcat(inputfilefolder,'Zmat_names'));
        for i = 1:N
            atom_names{i}  = zmat_names.textdata{i,3};
        end
    end

    
    %make atom names 2 characters long
    for i = 1:size(atom_names,1)
        if size(atom_names{i},2) == 1
            atom_names{i} = [atom_names{i}(1) ' '];
        end
    end
else
    atom_names = cell(1,size(atom_element, 1));
    for i = 1:size(atom_element, 1)
        atom_names{i} = horzcat(atom_element{i}, num2str(i));
    end
end

% Delete temporary files that have been produced
delete(horzcat(inputfilefolder,'Zmat_alone')); % Script makes this file
delete(horzcat(inputfilefolder,'Zmat_names')); % Script makes this file

end

