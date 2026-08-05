function [xyz_cg,xyz_Ca,chain_breakpoint,xyz_Ca_res,Ca]=pdbread_uaENM(filename_pdb)
[side_chain_type,res_pair_eps,res_pair_r0]=res_type_spica();
ft=fopen(filename_pdb,'r');
pdb = fread(ft);
fclose(ft);
lines_pdb = splitlines(string(char(pdb'))); 
pat="ATOM";
lines_pdb = lines_pdb(startsWith(lines_pdb,pat),:);
lines_pdb = char(lines_pdb);
xyz_cg = zeros(size(lines_pdb,1), 9);
xyz_cg(:,1:3) = [str2num(lines_pdb(:,31:38)), str2num(lines_pdb(:,39:46)), str2num(lines_pdb(:,47:54))];
chain=strip(string(lines_pdb(:,22)));
side_chain_Type=strip(string(lines_pdb(:,13:16)));
side_chain_Type(side_chain_Type=="TY1")="XYR";
side_chain_Type(side_chain_Type=="TY2")="BER";
side_chain_Type(side_chain_Type=="TY3")="PHR";
side_chain_Type(side_chain_Type=="TY4")="BER";
side_chain_Type(side_chain_Type=="PH1")="XYR";
side_chain_Type(side_chain_Type=="PH2")="BER";
side_chain_Type(side_chain_Type=="PH3")="BER";
side_chain_Type(side_chain_Type=="PH4")="BER";
side_chain_Type(side_chain_Type=="HS1")="HI1";
side_chain_Type(side_chain_Type=="HS2")="HI2";
side_chain_Type(side_chain_Type=="HS3")="HI3";
fesTyps=zeros(size(xyz_cg,1),1);
for i=1:size(side_chain_Type,1)  
     if isempty(find(matches(side_chain_type,side_chain_Type(i)), 1))
        disp(['There are unrecognizable residues',num2str(side_chain_Type(i)),'number',num2str(i)]);  
        error('Please check the input file.')
    end
    fesTyps(i)=find(matches(side_chain_type,side_chain_Type(i))); 
end
xyz_cg(:,4) = fesTyps;
%% 
Ca=[];
xyz_Ca=[];
xyz_Ca_res=[];
for i=1:size(lines_pdb,1)
    if side_chain_Type(i)=="GBM"
        xyz_Ca=[xyz_Ca;xyz_cg(i,1:3)];
        xyz_Ca_res= [xyz_Ca_res;strip(string(lines_pdb(i,17:20)))];
        Ca=[Ca;i];
    elseif side_chain_Type(i)=="ABB"
        xyz_Ca=[xyz_Ca;xyz_cg(i,1:3)];
        xyz_Ca_res= [xyz_Ca_res;strip(string(lines_pdb(i,17:20)))];
        Ca=[Ca;i];
    elseif side_chain_Type(i)=="GBB"
        xyz_Ca=[xyz_Ca;xyz_cg(i,1:3)];
        xyz_Ca_res= [xyz_Ca_res;strip(string(lines_pdb(i,17:20)))];
        Ca=[Ca;i];
    end
end

n_res=size(xyz_Ca,1);

res_number= [str2num(lines_pdb(:,23:26))];

res_number_new=zeros(size(res_number,1),1);
k_res=0;
for i=1:n_res-1  
  if side_chain_Type(k_res+1)=="ABB" || side_chain_Type(k_res+1)=="GBB"
       res_number_new(k_res+1)=i;
       k_res=k_res+1;
  else
    for j=1:6
        if res_number(k_res+j)==res_number(k_res+j+1)
            res_number_new(k_res+j)=i;
        elseif res_number(k_res+j)==res_number(k_res+j-1) 
            res_number_new(k_res+j)=i;
            k_res=k_res+j;
            break
        end
    end
  end
end

res_number_new(res_number_new==0)=n_res; 
chain_breakpoint=[];
for i=1:size(chain,1)-1
    if chain(i) ~= chain(i+1)
        chain_breakpoint=[chain_breakpoint;res_number_new(i)];    
    end
end
 
 xyz_cg(:,5)=res_number_new;

 elc=zeros(size(xyz_cg,1),1);
 elc(fesTyps==26)=0.559;
 elc(fesTyps==24)=0.559;
 elc(fesTyps==27)=-0.559;
 elc(fesTyps==28)=-0.559;
 xyz_cg(:,6) = elc;

 polar=zeros(size(xyz_cg,1),1);
 polar(fesTyps==19)=0.0559;
 polar(fesTyps==14)=0.0559;
 polar(fesTyps==20)=0.0559;
 polar(fesTyps==21)=0.0559;
 polar(fesTyps==22)=0.0559;
 polar(fesTyps==10)=0.0559;
 polar(fesTyps==13)=0.0559;
 xyz_cg(:,7) = polar;

vdw=diag(res_pair_eps);
vdw0=xyz_cg(:,4);
vdw0=vdw(vdw0);
xyz_cg(:,8)=vdw0;
r0=diag(res_pair_r0);
r00=xyz_cg(:,4);
r00=r0(r00);
xyz_cg(:,9)=r00;

