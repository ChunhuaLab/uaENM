clc
clear


mainDir = fileparts(mfilename('fullpath'));
parentDir = fileparts(mainDir);
addpath(mainDir, '-begin');
currentSubDir = fullfile(parentDir, 'example');
fileList_filename_pdb_cg = dir(fullfile(currentSubDir, '*-cg.pdb'));
assert(isscalar(fileList_filename_pdb_cg), 'uaENM:InputFiles', ...
    'Put exactly one protein in example: <id>-cg.pdb, <id>-y.pdb and <id>.top.');
proteinID = erase(fileList_filename_pdb_cg(1).name, '-cg.pdb');
fileNames_filename_pdb_cg = fullfile(currentSubDir, fileList_filename_pdb_cg(1).name);
fileNames_filename_pdb_o = fullfile(currentSubDir, [proteinID '-y.pdb']);
fileNames_filename_top = fullfile(currentSubDir, [proteinID '.top']);
assert(isfile(fileNames_filename_pdb_o), 'uaENM:MissingReference', ...
    'Missing experimental B-factor file: %s', fileNames_filename_pdb_o);
assert(isfile(fileNames_filename_top), 'uaENM:MissingTopology', ...
    'Missing topology file: %s', fileNames_filename_top);

[Smax_uaENM, cutoff, R_ab] = results_optimal_parameters( ...
    fileNames_filename_pdb_o, fileNames_filename_pdb_cg, fileNames_filename_top);
fprintf('The calculation results are as follows: pcc=%.12f, cutoff=%.2f, n_x=%.2f\n', ...
    Smax_uaENM, cutoff, R_ab);
resulPath = fullfile(currentSubDir, 'resul.txt');
fid = fopen(resulPath, 'w');
assert(fid >= 0, 'uaENM:OutputFile', 'Cannot write %s.', resulPath);
fprintf(fid, 'The calculation results are as follows: pcc=%.12f, cutoff=%.2f, n_x=%.2f\n', ...
    Smax_uaENM, cutoff, R_ab);
fclose(fid);
saveas(gcf, fullfile(currentSubDir, 'uaENM_result.png'));

function [Smax_uaENM, cutoff, R_ab] = results_optimal_parameters(filename_pdb_o, filename_pdb_cg, filename_top)

cutoffGrid = 10:20;
nrcGrid = 1.5:0.5:4.0;
num_R_ab = numel(nrcGrid);
num_jobs = numel(cutoffGrid)*num_R_ab;
all_S = NaN(num_jobs,1);
all_D_error = zeros(num_jobs,1);
all_warnings = zeros(num_jobs,1);
all_params = zeros(num_jobs,2);
[bond, angleparam, angle, dihedralparam, dist, xyz_cg, interaction, CA] = ...
    read_top(filename_pdb_cg, filename_top);
[~, res_pair_eps, res_pair_r0] = res_type_spica();
Smax_uaENM = -Inf;
cutoff = NaN;
R_ab = NaN;
bestResult = struct();

for idx = 1:num_jobs
    cutoff_temp = cutoffGrid(1 + floor((idx-1)/num_R_ab));
    R_ab_temp = nrcGrid(1 + mod(idx-1, num_R_ab));
    [S_CG_all, D_error, warning_cgall, n_temp, bf_temp, fitted_temp, msf_temp] = ...
        calculate_uaENM(filename_pdb_o, cutoff_temp, R_ab_temp, ...
        bond, angleparam, angle, dihedralparam, dist, xyz_cg, interaction, res_pair_eps, res_pair_r0, CA);
    all_S(idx) = S_CG_all(1,2);
    all_D_error(idx) = D_error;
    all_warnings(idx) = warning_cgall;
    all_params(idx,:) = [cutoff_temp, R_ab_temp];
    if all_D_error(idx) == 0 && all_warnings(idx) < 0.05 && isfinite(all_S(idx))
        if all_S(idx) > Smax_uaENM
            Smax_uaENM = all_S(idx);
            cutoff = all_params(idx,1);
            R_ab = all_params(idx,2);
            bestResult = struct('n',n_temp,'bf',bf_temp, ...
                'fitted',fitted_temp,'msf',msf_temp);
        end
    end
end
assert(isfinite(Smax_uaENM), 'uaENM:NoValidCombination', ...
    'No valid parameter combination: require lambda7 >= 1e-4 and capped fraction < 0.05.');


n = bestResult.n;
bf = bestResult.bf;
Bfactor_com = bestResult.fitted;
msf = bestResult.msf;
outputDir = fileparts(filename_pdb_o);
writematrix(msf, fullfile(outputDir, 'msf.txt'));
writematrix(Bfactor_com, fullfile(outputDir, 'Bfactor.txt'));
figure
aa = plot(1:n, bf, 'k', 'LineWidth', 1.2);
hold on
bb = plot(1:n, Bfactor_com, 'b', 'LineWidth', 1.2);
xlabel('Residue Index', 'FontSize', 15);
ylabel('B-factor (Å^2)', 'FontSize', 15);
legend([aa, bb], 'experimental', 'uaENM', 'FontSize', 13);
set(gcf, 'Color', 'white');
set(gca, 'LineWidth', 1.5, 'FontSize', 13, 'FontName', 'Times New Roman');
end
