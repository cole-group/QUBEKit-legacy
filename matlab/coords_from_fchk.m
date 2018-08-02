function [hessian,N] =  coords_from_fchk( inputfilefolder, fchk_file )
%Function extracts xyz file from the .fchk output file from Gaussian, this
%provides the coordinates of the molecules

if exist( horzcat(inputfilefolder,fchk_file), 'file') == 2
    fid = fopen(horzcat(inputfilefolder,fchk_file));
else
    fid_log = fopen(horzcat(outputfilefolder,'MSM_log'));
    fprintf(fid_log, '%s\n', 'ERROR - No .fchk file found.');
    fclose(fid_log);
end

tline = fgets(fid);

loop = 'y';
numbers = []; % Atomic numbers for use in xyz file
list_coords = []; % List of xyz coordinate
hessian =  [];

%Get atomic number and coordinates from fchk
while ischar(tline) && loop == 'y'
    %Gets atomic numbers
    if size(tline,2) > 14 && strcmp(strtrim(tline(1:14)), 'Atomic numbers')
        tline = fgets(fid);
        while size(tline,2) < 15 || strcmp(strtrim(tline(1:15)), 'Nuclear charges')  == 0
            tmp = strsplit(strtrim(tline));
            numbers(end + 1:(end + size(tmp,2))) = str2num(char(tmp));
            tline = fgets(fid);
        end
    end
    
    %Gets coordinates
    if size(tline,2) > 29 && strcmp(strtrim(tline(1:29)), 'Current cartesian coordinates')
        tline = fgets(fid);
        while size(tline,2) < 11 || strcmp(strtrim(tline(1:11)), 'Force Field')  == 0 && strcmp(strtrim(tline(1:14)), 'Int Atom Types')  == 0 && strcmp(strtrim(tline(1:10)), 'Atom Types')  == 0
            tmp = strsplit(strtrim(tline));
            list_coords(end + 1:(end + size(tmp,2))) = str2num(char(tmp)); 
            tline = fgets(fid);
        end
    end   
    
    %Gets Hessian
    if size(tline,2) > 25 && strcmp(strtrim(tline(1:25)), 'Cartesian Force Constants')
        tline = fgets(fid);
        while size(tline,2) < 13 || strcmp(strtrim(tline(1:13)), 'Dipole Moment')  == 0
            tmp = strsplit(strtrim(tline));
            hessian(end + 1:(end + size(tmp,2))) = str2num(char(tmp)); 
            tline = fgets(fid);
        end
        loop = 'n';
    end   
    
    tline = fgets(fid);
end


fclose(fid);

list_coords = list_coords.* 0.529; %Change from Bohr to Ang

N = size(list_coords,2)/3; %Number of atoms


%Opens .xyz file
fid = fopen(horzcat(inputfilefolder,'input_coords.xyz'), 'w'); 
fprintf(fid, '%1.0f \n \n',  N);
xyz = zeros(N,3);
n = 1;

names = cell(1, size(numbers,2));

%Change atomic numbers to atomic names
elements_numbers = importdata('elementlist.csv');
elements_names = cell(1, 99); 

for i =1:size(elements_numbers,1)
    elements_names{i} = elements_numbers{i}(4:5);
end

for i=1:size(numbers,2)
     names{i} = elements_names{numbers(i)};
end

%Print coordinates to new input_coords.xyz file
for i = 1:N
    for j = 1:3
        xyz(i,j) = list_coords(n);
        n = n + 1;
    end
    fprintf(fid, '%s  %4.4f   %4.4f  %4.4f \n',names{i}, xyz(i,1), xyz(i,2), xyz(i,3));

end

fclose(fid);

end

