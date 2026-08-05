clc
clear

mainDir = fileparts(mfilename('fullpath'));
parentDir = fileparts(mainDir);
addpath(genpath(parentDir));
currentSubDir = fullfile(parentDir, 'example');
fileList_filename_pdb_o = dir(fullfile(currentSubDir, '*-y.pdb'));
fileNames_filename_pdb_o = fileList_filename_pdb_o.name;
fileList_filename_pdb_cg = dir(fullfile(currentSubDir, '*-cg.pdb'));
fileNames_filename_pdb_cg = fileList_filename_pdb_cg.name;
fileList_filename_top = dir(fullfile(currentSubDir, '*.top'));
fileNames_filename_top = fileList_filename_top.name;

[Smax_uaENM,cutoff, R_ab]=results_optimal_parameters(fileNames_filename_pdb_o,fileNames_filename_pdb_cg,fileNames_filename_top);

fprintf('The calculation results are as follows：pcc=%.2f, cutoff=%.2f,n_x=%.2f \n',Smax_uaENM,cutoff, R_ab);

%% Save the fprintf output to resul.txt (located in the example directory)
resulPath = fullfile(currentSubDir, 'resul.txt');
fid = fopen(resulPath, 'w');
if fid > 0
    fprintf(fid, 'The calculation results are as follows：pcc=%.2f, cutoff=%.2f,n_x=%.2f \n', Smax_uaENM, cutoff, R_ab);
    fclose(fid);
end

%% Save the figure to the example directory
saveas(gcf, fullfile(currentSubDir, 'uaENM_result.png'));

%%
function [Smax_uaENM,cutoff, R_ab]=results_optimal_parameters(filename_pdb_o,filename_pdb_cg,filename_top)

Smax_uaENM = -10;
cutoff = 0;
R_ab = 0;
num_cutoffs = 14;
num_R_ab = 40;
num_jobs = num_cutoffs * num_R_ab;
all_S = zeros(num_jobs, 2); 
all_D_error = zeros(num_jobs, 1);
all_warnings = zeros(num_jobs, 1);
all_params = zeros(num_jobs, 2);
%% 
[bond, angleparam, angle, dihedralparam, dist, xyz_cg, interaction, CA] = ...
read_top(filename_pdb_cg, filename_top);
[~, res_pair_eps, res_pair_r0] = res_type_spica();
parfor idx = 1:num_jobs
    cutoff_temp = 5 + floor((idx-1) / num_R_ab);
    R_ab_temp = 0.1 * ((idx-1) - (cutoff_temp - 5) * num_R_ab + 1);
    [S_CG_all, D_error, warning_cgall,~,~,~] = ...
        calculate_uaENM(filename_pdb_o, cutoff_temp, R_ab_temp,...
        bond, angleparam, angle, dihedralparam, dist, ...
        xyz_cg, interaction, res_pair_eps, res_pair_r0, CA);
    all_S(idx, :) = [round(S_CG_all(1,2), 2), round(S_CG_all(1,1), 4)];
    all_D_error(idx) = D_error;
    all_warnings(idx) = warning_cgall;
    all_params(idx, :) = [cutoff_temp, R_ab_temp];
end
%% 
for idx = 1:num_jobs
    if all_D_error(idx) == 0 && all_warnings(idx) < 0.05
        current_S = all_S(idx, 1);
        if current_S > Smax_uaENM
            Smax_uaENM = current_S;
            cutoff = all_params(idx, 1);
            R_ab = all_params(idx, 2);
        end
    end
end

%% 
[~,~,~,n,bf,Bfactor_com] = calculate_uaENM(filename_pdb_o, cutoff, R_ab,...
         bond, angleparam, angle, dihedralparam, dist, ...
        xyz_cg, interaction, res_pair_eps, res_pair_r0, CA);
  
figure 
aa=plot(1:n,bf,'k','LineWidth',1.2);
hold on
bb=plot(1:n,Bfactor_com,'b','LineWidth',1.2);
hold on
xlabel('Residue Index','FontSize',15);
ylabel('Mean-square Fluctuation','FontSize',15);
legend([aa(1),bb(1)],'experimental','uaENM','FontSize',13);
set(gcf,'color','white');
set(gca,'linewidth',1.5); 
set(gca,'FontSize',13); 
set(gca,'FontName','Times New Rome'); 


end
