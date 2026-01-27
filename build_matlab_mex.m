% build_matlab_mex.m - Build MEX files for MATLAB interface
% Run this script from MATLAB in the verse directory to compile the MEX files

fprintf('Building VERSE MEX files...\n');

% Get the script directory
script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);

% Check if verse.c exists
if ~isfile('c/verse.c')
    error('verse.c not found in c/ subdirectory');
end

% Build MEX files
try
    fprintf('Compiling mintverse_mex.c...\n');
    mex('-R2017b', '-v', 'matlab/mintverse_mex.c', 'c/verse.c', '-o', 'matlab/mintverse_mex');
    fprintf('Successfully built mintverse MEX file\n\n');
catch ME
    if contains(ME.message, 'No supported compiler')
        fprintf('\nError: No C compiler configured.\n');
        fprintf('Run "mex -setup C" to configure a compiler, then try again.\n');
    else
        fprintf('Failed to build mintverse: %s\n', ME.message);
    end
    rethrow(ME);
end

try
    fprintf('Compiling minsarverse_mex.c...\n');
    mex('-R2017b', '-v', 'matlab/minsarverse_mex.c', 'c/verse.c', '-o', 'matlab/minsarverse_mex');
    fprintf('Successfully built minsarverse MEX file\n\n');
catch ME
    fprintf('Failed to build minsarverse: %s\n', ME.message);
    rethrow(ME);
end

fprintf('========================================\n');
fprintf('VERSE MATLAB interface built successfully!\n');
fprintf('Add matlab/ to your MATLAB path to use the functions.\n');
fprintf('========================================\n');
