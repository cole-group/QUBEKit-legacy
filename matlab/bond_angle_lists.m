function [bondlist,anglelist] =  bond_angle_lists(inputfilefolder)
%This function extracts a list of bond and angles from the Gaussian .log
%file

home = pwd;

cd(inputfilefolder)

%Finds the Gaussian .log file
if exist( 'zmat.log', 'file') == 2
    fid = fopen('zmat.log'); 
else
    if exist( 'lig.log', 'file') == 2
        fid = fopen('lig.log'); 
    else
        fid_log = fopen(horzcat(outputfilefolder,'MSM_log'));
        fprintf(fid_log, '%s\n', 'ERROR - No .log file found.');
        fclose(fid_log);
    end
end

tline = fgets(fid);
bondlist = zeros(1,2);
anglelist = zeros(1,2);

n = 1;
n_bond = 1;
n_angle = 1;
tmp = 'R';
B = [];

%Finds the bond and angles from the .log file
while ischar(tline)
    tline = fgets(fid);
    %Starts at point when bond and angle lists occur
    if size(tline,2) > 80 && strcmp(tline(1:81), ' ! Name  Definition              Value          Derivative Info.                !')
        tline = fgets(fid);
        tline = fgets(fid);
        %Stops when all bonds and angles recorded
        while ( strcmp( tmp(1), 'R') || strcmp( tmp(1), 'A') ) && strcmp(tline(2), '!') 
            
            A = strsplit(tline);
            tline = fgets(fid);
            
            B{n} = (A(4));
            tmp = char(B{n});
            
            %Bond List
            if  strcmp( tmp(1), 'R')
                tmp_splt = strsplit(tmp(3:(end - 1)), ',');
                bondlist(n_bond, 1) = str2double(tmp_splt{1});
                bondlist(n_bond, 2) = str2double(tmp_splt{2});
                n_bond = n_bond + 1;
            end
            
            %Angle List
            if  strcmp( tmp(1), 'A')
                tmp_splt = strsplit(tmp(3:(end - 1)), ',');
                anglelist(n_angle, 1) = str2double(tmp_splt{1});
                anglelist(n_angle, 2) = str2double(tmp_splt{2});
                anglelist(n_angle, 3) = str2double(tmp_splt{3});
                n_angle = n_angle + 1;
            end
            
            n = n + 1;
            
        end
        tline = -1; 
    end
end

cd(home)
 
end

