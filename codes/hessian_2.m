function [hessian2, warn_cg] = hessian_2(R_x, bond, twobody_pir, dist, xyz_cg)


% referenced calculation principle: Tirion, M. M.; Ben-Avraham, D. Phys. Rev. E 2015, 91, 32712.

siteCount = size(xyz_cg,1);
springLimit = 10;


if isempty(twobody_pir)
    twobody_pir = zeros(0,4);
end
pairRows = twobody_pir(:,1);
pairCols = twobody_pir(:,2);
wellDepth = twobody_pir(:,3);
referenceLength = twobody_pir(:,4);
pairDistance = dist(sub2ind(size(dist), pairRows, pairCols));


shiftedLength = ((13/7)^(1/6).*referenceLength)*R_x;
spring = ((72.*wellDepth).*(pairDistance-shiftedLength).^2) ...
    ./(pairDistance.^2.*(referenceLength-shiftedLength).^2);
warn_cg = sum(spring > springLimit)/numel(spring);
spring(spring < 0) = 0;
spring(spring > springLimit) = springLimit;

nonbonded = sparse(pairRows, pairCols, spring, siteCount, siteCount);
bonded = sparse(bond(:,1), bond(:,2), 50, siteCount, siteCount);
contactWeights = (bonded + bonded') + nonbonded;
hessian2 = assembleHessian(contactWeights, xyz_cg(:,1:3));
end

function hessian = assembleHessian(contactWeights, coordinates)

siteCount = size(contactWeights,1);
[siteRow, siteCol, weight] = find(tril(contactWeights));
separation = coordinates(siteCol,:) - coordinates(siteRow,:);
lengths = sqrt(sum(separation.^2,2));
direction = separation./lengths;


rowAxis = [1 1 1 2 2 2 3 3 3];
colAxis = [1 2 3 1 2 3 1 2 3];
projection = -direction(:,rowAxis).*direction(:,colAxis);
rowIndex = 3*siteRow - (3-rowAxis);
colIndex = 3*siteCol - (3-colAxis);
blockWeight = weight(:,ones(1,9));
hessian = sparse(rowIndex(:), colIndex(:), ...
    blockWeight(:).*projection(:), 3*siteCount, 3*siteCount);
hessian = hessian + hessian';


correction = blockSumDiag(hessian);
if ~isreal(hessian)

    hessian = hessian - correction;
    return
end
if isempty(hessian)
    hessian = hessian - full(correction);
    return
end
firstBlock = hessian(:,1:3) - full(correction(:,1:3));
assembled = zeros(3*siteCount,3*siteCount,'like',firstBlock);
assembled(:,1:3) = firstBlock;
for firstCol = 4:3:3*siteCount
    columns = firstCol:firstCol+2;
    assembled(:,columns) = hessian(:,columns) - full(correction(:,columns));
end
hessian = assembled;
end

function D = blockSumDiag(hessian)

axisSums = [sum(hessian(:,1:3:end),2), ...
    sum(hessian(:,2:3:end),2), sum(hessian(:,3:3:end),2)];
cartesianSize = size(axisSums,1);
if ~isreal(hessian)

    D = zeros(cartesianSize);
    for firstRow = 1:3:cartesianSize
        blockRows = firstRow:firstRow+2;
        D(blockRows,blockRows) = axisSums(blockRows,:);
    end
    return
end

rowIndex = repmat((1:cartesianSize)',1,3);
colIndex = 3*floor(((1:cartesianSize)'-1)/3) + (1:3);
D = sparse(rowIndex(:),colIndex(:),double(axisSums(:)),cartesianSize,cartesianSize);
end

