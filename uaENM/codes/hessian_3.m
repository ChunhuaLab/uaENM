function [hessian3]=hessian_3(angleparam,angle,xyz_cg)
angles=[angleparam;angle];
%% 
ktheta = 20;
ktheta = ktheta(ones(size(angles,1),1)); 
%% 
x=xyz_cg(:,1:3);
clear xyz_cg;
%%
hessian3=zeros(size(x,1)*3);
for m=1:size(angles,1)
	K_theta = ktheta(m); 
    r=x(angles(m,:),1:3);
    A = r(1,:)-r(2,:);
	B = r(3,:)-r(2,:);
	G = cross(A,B);
	G = G/norm(G);
	dphidAB = [cross(G,A)/sum(A.^2); cross(B,G)/sum(B.^2);]'; 
	dABdr = [1 -1 0; 
			 0 -1 1;]; 
	dphiDr = dphidAB*dABdr;
	dphi2 = dphiDr(:)*dphiDr(:)';
	
	idx = [3*(angles(m,:)-1)+1; 3*(angles(m,:)-1)+2; 3*(angles(m,:)-1)+3];
    hessian3(idx,idx) = hessian3(idx,idx) + (2*K_theta)*dphi2; 
end
hessian3=sparse(hessian3);
