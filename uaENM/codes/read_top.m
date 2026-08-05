function [bond,angleparam,angle,dihedralparam,dist,xyz_cg,interaction,CA]=read_top(filename_pdb,filename_top)
% Read topology file and build the non-bonded interaction matrix.
% Note: 'bondparam' lines are force-field parameters (not bond records) and are
% skipped, so they are not misread as 'bond' entries by startsWith(line,'bond').

[xyz_cg,~,~,~,CA]=pdbread_uaENM(filename_pdb);
dist = squareform(pdist(xyz_cg(:,1:3)));
bond = [];
angleparam = [];
angle = [];
dihedralparam = [];
%% 
    fid = fopen(filename_top, 'r');
    while ~feof(fid)
        line = strtrim(fgetl(fid)); 
        if isempty(line)
            continue; 
        end
        if startsWith(line, 'bondparam')
            continue; % bondparam lines are force-field parameters, not bonds; skip
        elseif startsWith(line, 'bond')
            parts = split(line);
            nums = str2double(parts(2:3)); 
            if all(~isnan(nums))
                bond = [bond; nums'];
            end
        elseif startsWith(line, 'angleparam')
            parts = split(line);
            nums = str2double(parts(2:4)); 
            if all(~isnan(nums))
                angleparam = [angleparam; nums'];
            end
        elseif startsWith(line, 'angle') && ~contains(line, 'angleparam')
            parts = split(line);
            nums = str2double(parts(2:4)); 
            if all(~isnan(nums))
                angle = [angle; nums'];
            end
        elseif startsWith(line, 'dihedralparam')
            parts = split(line);
            nums = str2double(parts(2:5)); 
            if all(~isnan(nums))
                dihedralparam = [dihedralparam; nums'];
            end
        end
    end
    
    fclose(fid);
%% 
   interaction = zeros(size(dist,1));
  for i = 1:size(bond, 1)
    a1 = bond(i, 1);
    a2 = bond(i, 2);
    interaction(a1, a2) = 2;
    interaction(a2, a1) = 2;
  end    

  angles=[angleparam;angle];
  for i = 1:size(angles, 1)
    a1 = angles(i, 1);
    a2 = angles(i, 2);
    a3 = angles(i, 3);
    interaction(a1, a2) = 3;
    interaction(a2, a1) = 3;
    interaction(a2, a3) = 3;
    interaction(a3, a2) = 3;
    interaction(a1, a3) = 3;
    interaction(a3, a1) = 3;
  end
 for i = 1:size(dihedralparam, 1)
    a1 = dihedralparam(i, 1);
    a2 = dihedralparam(i, 2);
    a3 = dihedralparam(i, 3);
    a4 = dihedralparam(i, 4);
    pairs = [a1 a2; a1 a3; a1 a4; a2 a3; a2 a4;a3 a4];
    for j = 1:size(pairs, 1)
        p1 = pairs(j, 1);
        p2 = pairs(j, 2);
        interaction(p1, p2) = 4;
        interaction(p2, p1) = 4;
    end
 end 
end
