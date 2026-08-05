function [S,D_error,warn_cg,n,bf,Bfactor_com]=calculate_uaENM(filename_pdb_bf,cutoff,R_x,...
 bond,angleparam,angle,dihedralparam,dist,xyz_cg,interaction,res_pair_eps,res_pair_r0 ,CA)

[twobody_pir]=confirm_twobody_pir(res_pair_eps,res_pair_r0,interaction,dist,xyz_cg,cutoff);
[hessian2,warn_cg]=hessian_2(R_x,bond,twobody_pir,dist,xyz_cg);
hessian3=hessian_3(angleparam,angle,xyz_cg);
hessian4=hessian_4(dihedralparam,xyz_cg);
hessian=hessian2+hessian3+hessian4;
hessian=full(hessian);

[V,D]=eig(hessian);
clear hessian

if D(7,7) < 1e-4
    D_error =1;
else
    D_error = 0;
end
n = size(CA, 1);

row_indices = reshape(3*CA' - [2; 1; 0], 1, [])';
V_CA = V(row_indices, :);

rows1 = 3*(1:n) - 2; 
rows2 = 3*(1:n) - 1;
rows3 = 3*(1:n);
cols_idx = 7:size(V_CA,2);  
V1 = V_CA(rows1, cols_idx);    
V2 = V_CA(rows2, cols_idx);      
V3 = V_CA(rows3, cols_idx);    
D_diag = diag(D(cols_idx, cols_idx))';  
msf = sum((V1.^2 + V2.^2 + V3.^2) ./ D_diag, 2); 


ft=fopen(filename_pdb_bf,'r');
pdb = fread(ft);
fclose(ft);
lines_pdb = splitlines(string(char(pdb')));
pat="ATOM";
lines_pdb = lines_pdb(startsWith(lines_pdb,pat),:);
lines_pdb = char(lines_pdb);
atom_type=strip(string(lines_pdb(:,13:16)));
bf_temp = [str2num(lines_pdb(:,61:66))];
bf=[];
for i=1:size(lines_pdb,1)
    if atom_type(i)=="CA"
        bf=[bf;bf_temp(i)];
    end
end

cor(:,2)=msf;
cor(:,1)=bf;

[S,P]=corr(cor,'type','Pearson');
S(1,1)=P(1,2);
S(2,2)=P(2,1);
warnState = warning('off', 'stats:regress:RankDefDesignMat');
B=regress(bf,[ones(n,1) msf]);
warning(warnState);
for i = 1:n
    Bfactor_com(i) = (B(1)+B(2)*msf(i))';
end
