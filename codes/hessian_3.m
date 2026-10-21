function hessian3 = hessian_3(angleparam, angle, xyz_cg)


% referenced calculation principle: Blondel and Karplus, J. Comput. Chem., Vol. 17, No. 9, 1132-1141 (1996)
% reference code: Na, H.; Song, G., Phys. Biol. 2023, 20, 46005.

triplets = [angleparam; angle];
coordinates = xyz_cg(:,1:3);
angleStiffness = 20;
edgeToSite = [1 -1 0; 0 -1 1];
cartesianSize = 3*size(coordinates,1);
if isreal(coordinates)
    hessian3 = spalloc(cartesianSize,cartesianSize, ...
        min(cartesianSize^2,81*size(triplets,1)));
else
    hessian3 = zeros(cartesianSize);
end


for record = 1:size(triplets,1)
    sites = triplets(record,:);
    points = coordinates(sites,1:3);
    leftEdge = points(1,:) - points(2,:);
    rightEdge = points(3,:) - points(2,:);
    planeNormal = cross(leftEdge, rightEdge);
    planeNormal = planeNormal/norm(planeNormal);

    leftDerivative = cross(planeNormal, leftEdge)/sum(leftEdge.^2);
    rightDerivative = cross(rightEdge, planeNormal)/sum(rightEdge.^2);
    edgeGradient = [leftDerivative; rightDerivative]';
    siteGradient = edgeGradient*edgeToSite;
    curvature = siteGradient(:)*siteGradient(:)';

    cartesian = [3*(sites-1)+1; 3*(sites-1)+2; 3*(sites-1)+3];

    previousBlock = full(hessian3(cartesian,cartesian));

    hessian3(cartesian,cartesian) = previousBlock ...
        + (2*angleStiffness)*curvature;
end
hessian3 = sparse(hessian3);
end
