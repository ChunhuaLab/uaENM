function [twobody_pir]=confirm_twobody_pir(res_pair_eps,res_pair_r0,interaction,dist,xyz_cg,cutoff)

[rows, cols] = find(interaction == 0);
pairs = [rows, cols];
dists = dist(sub2ind(size(dist), rows, cols)); 
valid = (rows ~= cols) & (dists < cutoff);   
twobody_pir_0 = pairs(valid, :);
%% 
if ~isempty(twobody_pir_0)
    type1 = xyz_cg(twobody_pir_0(:,1), 4);
    type2 = xyz_cg(twobody_pir_0(:,2), 4);
    indices = sub2ind(size(res_pair_eps), type1, type2);
    eps_values = res_pair_eps(indices);
    r0_values = res_pair_r0(indices);
    twobody_pir = [twobody_pir_0, eps_values, r0_values];
else
    twobody_pir = [];
end
end