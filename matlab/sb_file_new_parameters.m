function  sb_file_new_parameters(inputfilefolder, filename)
%Takes new angle and bond terms and puts them into a .sb file with name
%filename_seminario.sb

delimiterIn = ',';

angles  = importdata(horzcat(inputfilefolder,'Average_Modified_Seminario_Angles'), delimiterIn); 

bonds  = importdata(horzcat(inputfilefolder,'Average_Modified_Seminario_Bonds'), delimiterIn);

reversed_angles = angles; 

%Formats the angle terms
for i =1:size(angles.textdata,1)
    tmp = angles.textdata{i}; 
    %Makes 8 characters long 
    if size(tmp) < 8 
        tmp = horzcat(tmp, ' '); 
    end
    reversed_angles.textdata{i} = horzcat(tmp(7:8), tmp(3:6), tmp(1:2)); 
end

%Opens Files
fidout = fopen(horzcat(inputfilefolder,filename,'_Seminario.sb'), 'wt'); % Script makes this file
fprintf(fidout, '*****                         Bond Stretching and Angle Bending Parameters - July 17*****\n');

%Writes to file

%Prints out bonds at top of file
for i = 1:(size(bonds.textdata,1))
    if size(bonds.textdata{i}, 2) == 5
        fprintf(fidout, '%s %3.1f      %3.3f        Modified Seminario Method AEAA \n', bonds.textdata{i}, bonds.data(i,1), bonds.data(i,2));
    else
        fprintf(fidout, '%s  %3.1f      %3.3f        Modified Seminario Method AEAA \n', bonds.textdata{i}, bonds.data(i,1), bonds.data(i,2));
    end
end

fprintf(fidout, '\n');
fprintf(fidout, '********                        line above must be blank\n');

%Prints out angles in middle section of file
for i = 1:(size(angles.textdata,1))
    if size(angles.textdata{i},2) == 8
        fprintf(fidout, '%s    %3.1f       %3.1f    Modified Seminario Method AEAA \n', angles.textdata{i}, angles.data(i,1), angles.data(i,2));
    else
        fprintf(fidout, '%s     %3.1f       %3.1f    Modified Seminario Method AEAA \n', angles.textdata{i}, angles.data(i,1), angles.data(i,2));
    end
end

fprintf(fidout, '\n');
fprintf(fidout, '\n');

fclose(fidout);
end