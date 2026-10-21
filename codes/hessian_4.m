function hessian4 = hessian_4(phis, xyz_cg)

% referenced calculation principle: Blondel and Karplus, J. Comput. Chem., Vol. 17, No. 9, 1132-1141 (1996)
% reference code: Na, H.; Song, G., Phys. Biol. 2023, 20, 46005.

coordinates = xyz_cg(:,1:3);
if isempty(phis)

    hessian4 = zeros(3*size(coordinates,1), 3*size(coordinates,1));
else
    recordCount = size(phis,1);
    amplitude = 50*ones(recordCount,1);
    periodicity = ones(recordCount,1);
    parameterized = [phis, amplitude, periodicity];
    stiffness = parameterized(:,5).*(parameterized(:,6)).^2/2;
    hessian4 = STeMphi(coordinates, parameterized(:,1:4), stiffness);
end

    function hessians = STeMphi(coordinates, quadruplets, stiffness)

        edgeToSite = [1 -1 0 0; 0 1 -1 0; 0 0 -1 1];
        cartesianSize = 3*size(coordinates,1);
        if isreal(coordinates)
            accumulator = spalloc(cartesianSize,cartesianSize, ...
                min(cartesianSize^2,144*size(quadruplets,1)));
        else
            accumulator = zeros(cartesianSize);
        end
        for record = 1:size(quadruplets,1)
            sites = quadruplets(record,:);
            points = coordinates(sites,1:3);
            firstEdge = points(1,:) - points(2,:);
            middleEdge = points(2,:) - points(3,:);
            lastEdge = points(4,:) - points(3,:);
            firstNormal = cross(firstEdge, middleEdge);
            lastNormal = cross(lastEdge, middleEdge);
            firstArea2 = sum(firstNormal.^2);
            lastArea2 = sum(lastNormal.^2);
            middleLength = norm(middleEdge);

            firstDerivative = -middleLength/firstArea2*firstNormal;
            middleDerivative = dot(firstEdge,middleEdge)/firstArea2/middleLength*firstNormal ...
                - dot(lastEdge,middleEdge)/lastArea2/middleLength*lastNormal;
            lastDerivative = middleLength/lastArea2*lastNormal;
            edgeGradient = [firstDerivative; middleDerivative; lastDerivative]';
            siteGradient = edgeGradient*edgeToSite;
            curvature = siteGradient(:)*siteGradient(:)';

            cartesian = [3*(sites-1)+1; 3*(sites-1)+2; 3*(sites-1)+3];

            previousBlock = full(accumulator(cartesian,cartesian));

            accumulator(cartesian,cartesian) = previousBlock ...
                + (2*stiffness(record))*curvature;
        end
        hessians = sparse(accumulator);
    end
end
