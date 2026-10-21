function [S, D_error, warn_cg, n, bf, Bfactor_com, msf] = calculate_uaENM(filename_pdb_bf, cutoff, R_x, ...
    bond, angleparam, angle, dihedralparam, dist, xyz_cg, interaction, res_pair_eps, res_pair_r0, CA)


persistent geometryCache referenceCache
twobody_pir = confirm_twobody_pir(res_pair_eps, res_pair_r0, interaction, dist, xyz_cg, cutoff);
[hessian2, warn_cg] = hessian_2(R_x, bond, twobody_pir, dist, xyz_cg);


geometryInputs = {xyz_cg(:,1:3), angleparam, angle, dihedralparam};
if isempty(geometryCache) || ~sameGeometryInputs(geometryInputs,geometryCache.inputs)
    hessian3 = hessian_3(angleparam, angle, xyz_cg);
    hessian4 = hessian_4(dihedralparam, xyz_cg);
    geometryCache = struct('inputs',{geometryInputs},'hessian3',hessian3,'hessian4',hessian4);
else
    hessian3 = geometryCache.hessian3;
    hessian4 = geometryCache.hessian4;
end
hessian = full(hessian2 + hessian3 + hessian4);
assert(all(isfinite(hessian(:))), 'uaENM:InvalidHessian', ...
    'Nonfinite Hessian entries; check coincident atoms and degenerate geometry.');
[V, D] = eig(hessian);
[D_diag, order] = sort(diag(D));
V = V(:,order);
assert(numel(D_diag) >= 7, 'uaENM:TooFewAtoms', 'At least seven Cartesian modes are required.');
D_error = double(D_diag(7) < 1e-4);
n = size(CA,1);
row_indices = reshape(3*CA' - [2; 1; 0], 1, [])';
V_CA = V(row_indices,:);
cols_idx = 7:size(V_CA,2);
V1 = V_CA(3*(1:n)-2, cols_idx);
V2 = V_CA(3*(1:n)-1, cols_idx);
V3 = V_CA(3*(1:n), cols_idx);
msf = sum((V1.^2 + V2.^2 + V3.^2) ./ D_diag(cols_idx)', 2);


referenceText = fileread(filename_pdb_bf);
if isempty(referenceCache) || ~isequal(referenceText,referenceCache.text)
    lines = splitlines(string(referenceText));
    lines = lines(startsWith(lines, 'ATOM'));
    assert(~isempty(lines) && all(strlength(lines) >= 66), 'uaENM:InvalidReference', ...
        'Experimental PDB must contain complete ATOM records including B-factor columns 61-66.');
    atoms = char(lines);
    bf = str2double(string(atoms(strtrim(string(atoms(:,13:16))) == "CA", 61:66)));
    referenceCache = struct('text',referenceText,'bf',bf);
else
    bf = referenceCache.bf;
end
assert(numel(bf) == n && all(isfinite(bf)) && norm(bf-mean(bf)) > 0, ...
    'uaENM:InvalidReference', 'Reference CA count must match the CG structure; B-factors must be finite and nonconstant.');

S = NaN(2);
if all(isfinite(msf)) && norm(msf-mean(msf)) > 0
    x = bf - mean(bf);
    y = msf - mean(msf);
    pcc = (x'*y)/(norm(x)*norm(y));
    S(1,2) = pcc;
    S(2,1) = pcc;
    B = [ones(n,1), msf] \ bf;
    Bfactor_com = B(1) + B(2)*msf;
else
    Bfactor_com = NaN(n,1);
end
end

function same = sameGeometryInputs(current,previous)

same = isequaln(current,previous);
if ~same, return; end
for k = 1:numel(current)
    a = current{k}; b = previous{k};
    if ~strcmp(class(a),class(b)) || issparse(a)~=issparse(b) || isreal(a)~=isreal(b)
        same = false;
        return
    end
    if isa(a,'double') || isa(a,'single')
        a = full(a(:)); b = full(b(:));
        if ~isequal(typecast(real(a),'uint8'),typecast(real(b),'uint8')) ...
                || ~isequal(typecast(imag(a),'uint8'),typecast(imag(b),'uint8'))
            same = false;
            return
        end
    end
end
end
