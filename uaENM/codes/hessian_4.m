function [hessian4]=hessian_4(phis,xyz_cg)

% referenced calculation principle: Blondel and Karplus, J. Comput. Chem., Vol. 17, No. 9, 1132-1141 (1996)
% reference code: Na, H.; Song, G., Phys. Biol. 2023, 20, 46005.

x=xyz_cg(:,1:3);
if isempty(phis)
    hessian4=zeros(size(x,1)*3, size(x,1)*3);
else
Kphi_K=50;
Kphi_K = Kphi_K(ones(size(phis,1),1)); 
Kphi_n = 1;
Kphi_n = Kphi_n(ones(size(phis,1),1)); 
Phis_K = [Kphi_K,Kphi_n];
kphi = [phis, Phis_K];
hessian4=STeMphi(x,kphi(:,1:4),kphi(:,5).*(kphi(:,6)).^2/2);
end

function [hessians] =STeMphi(x,torsional,kphi)
hessian=zeros(size(x,1)*3, size(x,1)*3);
for m=1:size(torsional,1)
	K_phi = kphi(m); 
    r=x(torsional(m,:),1:3);
    F = r(1,:)-r(2,:);
	G = r(2,:)-r(3,:);
	H = r(4,:)-r(3,:);
	A = cross(F,G);
	B = cross(H,G);
	A2 = sum(A.^2);
	B2 = sum(B.^2);
	dphidFGH = [ -norm(G)/A2*A;	
			dot(F,G)/A2/norm(G)*A - dot(H,G)/B2/norm(G)*B;
			norm(G)/B2*B]'; 
	dFGHdr = [1 -1 0 0; 
			  0 1 -1 0; 
			  0 0 -1 1]; 
   	dphiDr = dphidFGH*dFGHdr;
	dphi2 = dphiDr(:)*dphiDr(:)';
	idx_ST = [3*(torsional(m,:)-1)+1; 3*(torsional(m,:)-1)+2; 3*(torsional(m,:)-1)+3];
	hessian(idx_ST,idx_ST) = hessian(idx_ST,idx_ST) + (2*K_phi)*dphi2; 
end
     hessians=sparse(hessian);
end

end
