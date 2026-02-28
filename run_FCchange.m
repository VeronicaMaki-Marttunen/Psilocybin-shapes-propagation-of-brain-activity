% Calculates FC distance  matrix

numRegions=200;

cont1 = 1;
cont2 = 1;
cont3 = 1;
cont4 = 1;
cont5 = 1;
for i =1:size(connAll,1)
    Mx_z_sparse = squeeze(connAll(i,1,:,:));
    symmetric_matrix = Mx_z_sparse + Mx_z_sparse';
    symmetric_matrix(eye(size(symmetric_matrix)) == 1) = 0;
    non_zero_indices = any(symmetric_matrix, 2);
    filtered_matrix = symmetric_matrix(non_zero_indices, non_zero_indices);
    if ~isempty(filtered_matrix)
        connmatrixAll_B1{i} = filtered_matrix;
        cont1 = cont1 + 1;
    end
    clear normalized_matrix
    Mx_z_sparse = squeeze(connAll(i,2,:,:));
    symmetric_matrix = Mx_z_sparse + Mx_z_sparse';
    symmetric_matrix(eye(size(symmetric_matrix)) == 1) = 0;
    non_zero_indices = any(symmetric_matrix, 2);
    filtered_matrix = symmetric_matrix(non_zero_indices, non_zero_indices);
    if ~isempty(filtered_matrix)
        connmatrixAll_B2{i} = filtered_matrix;
        cont2 = cont2 + 1;
    end
    clear normalized_matrix

    Mx_z_sparse = squeeze(connAll(i,3,:,:));
    symmetric_matrix = Mx_z_sparse + Mx_z_sparse';
    symmetric_matrix(eye(size(symmetric_matrix)) == 1) = 0;
    non_zero_indices = any(symmetric_matrix, 2);
    filtered_matrix = symmetric_matrix(non_zero_indices, non_zero_indices);
    if ~isempty(filtered_matrix)
        connmatrixAll_B3{i} = filtered_matrix;
        cont3 = cont3 + 1;
    end
    clear normalized_matrix

    Mx_z_sparse = squeeze(connAll(i,4,:,:));
    symmetric_matrix = Mx_z_sparse + Mx_z_sparse';
    symmetric_matrix(eye(size(symmetric_matrix)) == 1) = 0;
    non_zero_indices = any(symmetric_matrix, 2);
    filtered_matrix = symmetric_matrix(non_zero_indices, non_zero_indices);
    if ~isempty(filtered_matrix)

        connmatrixAll_D{i} = filtered_matrix;
        cont4 = cont4 + 1;;
    end
    clear normalized_matrix

    Mx_z_sparse = squeeze(connAll(i,5,:,:));
    symmetric_matrix = Mx_z_sparse + Mx_z_sparse';
    symmetric_matrix(eye(size(symmetric_matrix)) == 1) = 0;
    non_zero_indices = any(symmetric_matrix, 2);
    filtered_matrix = symmetric_matrix(non_zero_indices, non_zero_indices);
    if ~isempty(filtered_matrix)

        connmatrixAll_M{i} = filtered_matrix;
        cont5 = cont5 + 1;;
    end
    clear normalized_matrix


end

% Extract upper triangles
for i = 1:7
    clear distance
    mx1 = connmatrixAll_B1{i};
    mx2 = connmatrixAll_B2{i};
    mx3 = connmatrixAll_B3{i};
    mx4 = connmatrixAll_D{i};
    mx5 = connmatrixAll_M{i};
    if length(mx3>1)
        upperTri1 = mx1(triu(true(numRegions), 1)); % Indices of upper triangular part
        upperTri2 = mx2(triu(true(numRegions), 1));
        upperTri3 = mx3(triu(true(numRegions), 1));
        upperTri4 = mx4(triu(true(numRegions), 1));
        upperTri5 = mx5(triu(true(numRegions), 1));
        % Calculate the RMS Euclidean distance
        distance(1,1) = sqrt(mean((upperTri1 - upperTri2).^2));
        distance(1,2) = sqrt(mean((upperTri1 - upperTri3).^2));
        distance(1,3) = sqrt(mean((upperTri1 - upperTri5).^2));
        distance(1,4) = sqrt(mean((upperTri1 - upperTri4).^2));
        distance(2,2) = sqrt(mean((upperTri2 - upperTri3).^2));
        distance(2,3) = sqrt(mean((upperTri2 - upperTri5).^2));
        distance(2,4) = sqrt(mean((upperTri2 - upperTri4).^2));
        distance(3,3) = sqrt(mean((upperTri3 - upperTri5).^2));
        distance(3,4) = sqrt(mean((upperTri3 - upperTri4).^2));
        distance(4,4) = sqrt(mean((upperTri5 - upperTri4).^2));
    else
        upperTri1 = mx1(triu(true(numRegions), 1));
        upperTri2 = mx2(triu(true(numRegions), 1));
        upperTri4 = mx4(triu(true(numRegions), 1));
        distance(1,1) = sqrt(mean((upperTri1 - upperTri2).^2));
        distance(1,3) = sqrt(mean((upperTri1 - upperTri4).^2));
        distance(2,3) = sqrt(mean((upperTri2 - upperTri4).^2));
        distance(3,3) = NaN;
        distance(3,4) = NaN;
    end
    % Store the distance in the distance matrix
    distanceMatrix{i} = distance;
end

for i = 1:7
    distce(i,:,:) = distanceMatrix{i};
end
figure
imagesc(squeeze(nanmean(distce)))
set(gca,'ytick',1:4,'yticklabel',{'B1','B2','B3','M'})
set(gca,'xtick',1:4,'xticklabel',{'B2','B3','M','PSI'})
set(gca,'fontsize',14)
colormap('gray')