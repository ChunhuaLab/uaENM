function [hessian2,warn_cg]=hessian_2(R_x,bond,twobody_pir,dist,xyz_cg)

% referenced calculation principle: Blondel and Karplus, J. Comput. Chem., Vol. 17, No. 9, 1132-1141 (1996)
% reference code: Na, H.; Song, G., Phys. Biol. 2023, 20, 46005.

%% vdW
n = size(xyz_cg,1);
eps=twobody_pir(:,3);
r0=twobody_pir(:,4);
r = dist(sub2ind(size(dist), twobody_pir(:,1), twobody_pir(:,2)));
R_ab=((13/7)^(1/6).*r0)*R_x;
k_vdw=((72.*eps).*(r-R_ab).^2)./(r.^2.*(r0-R_ab).^2);
%elc
k_elc=664*(xyz_cg(twobody_pir(:,1),6).*xyz_cg(twobody_pir(:,2),6))./r.^3;
k_all=k_vdw+k_elc;
warn_cg=size(find(k_all>14.29),1)/size(k_all,1);
k_all(k_all<0)=0;
k_all(k_all>14.29)=14.29;
cx_vdW = sparse(twobody_pir(:,1), twobody_pir(:,2), k_all, n, n);
%% bond
cx_bond = sparse(bond(:,1),bond(:,2),50,n,n);
clear bds bonds;
%% 
cx_sum=(cx_bond + cx_bond') + cx_vdW ;
xyz_r=xyz_cg(:,1:3);
cx_sum = tril(cx_sum);
dim = size(cx_sum,1);
[I,J,K] = find(cx_sum);
dR  = xyz_r(J,:) - xyz_r(I,:);
dr = sqrt(sum(dR.^2,2));
clear cx_UB cx_bond cx_sum cx_vdW  xyz xyz_r;
dR = dR./repmat(dr, 1, 3);
dRa = [repmat(dR(:,1),1,3), repmat(dR(:,2),1,3), repmat(dR(:,3),1,3)]; 
dRb = repmat(dR, 1, 3);
dR2 = -dRa.*dRb; 
Is = [repmat(3*I-2, 1, 3), repmat(3*I-1, 1, 3), repmat(3*I, 1, 3)];
Js = repmat([3*J-2, 3*J-1, 3*J], 1, 3);
Ks = repmat(K, 1, 9);
hessian2 = sparse(Is(:), Js(:), Ks(:).*dR2(:), 3*dim, 3*dim);
hessian2 = hessian2 + hessian2';
hessian2 = hessian2 - blockSumDiag(hessian2);
%% 
function [D] = blockSumDiag(hv1)
blockSum = [sum(hv1(:,1:3:end),2), sum(hv1(:,2:3:end),2), sum(hv1(:,3:3:end),2)];
n_b = size(blockSum,1);
D = zeros(n_b);
for i_b=1:n_b/3
	D((i_b-1)*3+1:i_b*3, (i_b-1)*3+1:i_b*3) = blockSum((i_b-1)*3+1:i_b*3,1:3);
end
end
end
