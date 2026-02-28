% addpath BrainSpace toolbox to path
% Schaefer parcellation file should be on path

load conn_Schaef
[surf_lh, surf_rh] = load_conte69();


parcelSch=ft_read_cifti('Schaefer2018_200Parcels_7Networks_order.dlabel.nii');
labeling = parcelSch.parcels(1:64984);

% All_subjects
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
    % Normalize each row of the filtered matrix
    if ~isempty(filtered_matrix)
        norms = sqrt(sum(filtered_matrix.^2, 2));
        normalized_matrix = filtered_matrix ./ norms;

        connmatrixAll_B1{cont1} = normalized_matrix;
        cont1 = cont1 + 1;
    end
    clear normalized_matrix
    Mx_z_sparse = squeeze(connAll(i,2,:,:));
    symmetric_matrix = Mx_z_sparse + Mx_z_sparse';
    symmetric_matrix(eye(size(symmetric_matrix)) == 1) = 0;
    non_zero_indices = any(symmetric_matrix, 2);
    filtered_matrix = symmetric_matrix(non_zero_indices, non_zero_indices);
    % Normalize each row of the filtered matrix
    if ~isempty(filtered_matrix)
        norms = sqrt(sum(filtered_matrix.^2, 2));
        normalized_matrix = filtered_matrix ./ norms;

        connmatrixAll_B2{cont2} = normalized_matrix;
        cont2 = cont2 + 1;
    end
    clear normalized_matrix

    Mx_z_sparse = squeeze(connAll(i,3,:,:));
    symmetric_matrix = Mx_z_sparse + Mx_z_sparse';
    symmetric_matrix(eye(size(symmetric_matrix)) == 1) = 0;
    non_zero_indices = any(symmetric_matrix, 2);
    filtered_matrix = symmetric_matrix(non_zero_indices, non_zero_indices);
    if ~isempty(filtered_matrix) % Check if there are non-zero rows/columns
        norms = sqrt(sum(filtered_matrix.^2, 2));
        normalized_matrix = filtered_matrix ./ norms;
        connmatrixAll_B3{cont3} = normalized_matrix;
        cont3 = cont3 + 1;
    end
    clear normalized_matrix

    Mx_z_sparse = squeeze(connAll(i,4,:,:));
    symmetric_matrix = Mx_z_sparse + Mx_z_sparse';
    symmetric_matrix(eye(size(symmetric_matrix)) == 1) = 0;
    non_zero_indices = any(symmetric_matrix, 2);
    filtered_matrix = symmetric_matrix(non_zero_indices, non_zero_indices);
    if ~isempty(filtered_matrix) % Check if there are non-zero rows/columns
        norms = sqrt(sum(filtered_matrix.^2, 2));
        normalized_matrix = filtered_matrix ./ norms;

        connmatrixAll_D{cont4} = normalized_matrix;
        cont4 = cont4 + 1;
    end
    clear normalized_matrix


    Mx_z_sparse = squeeze(connAll(i,5,:,:));
    symmetric_matrix = Mx_z_sparse + Mx_z_sparse';
    symmetric_matrix(eye(size(symmetric_matrix)) == 1) = 0;
    non_zero_indices = any(symmetric_matrix, 2);
    filtered_matrix = symmetric_matrix(non_zero_indices, non_zero_indices);
    if ~isempty(filtered_matrix) % Check if there are non-zero rows/columns
        norms = sqrt(sum(filtered_matrix.^2, 2));
        normalized_matrix = filtered_matrix ./ norms;

        connmatrixAll_M{cont5} = normalized_matrix;
        cont5 = cont5 + 1;;
    end
    clear normalized_matrix

end
Gp = GradientMaps('kernel','cs','approach','diffusionEmbedding','alignment','pa','n_components',5);

% Calculate gradient scores

Gp = Gp.fit(connmatrixAll_B1);
for i= 1:11
    pg0(i,:)=Gp.aligned{i}(:,1);
end

Gp = Gp.fit(connmatrixAll_B2);
for i= 1:7
    pg(i,:)=Gp.aligned{i}(:,1);
end

Gp = Gp.fit(connmatrixAll_B3);
for i= 1:7
    pg2(i,:)=Gp.aligned{i}(:,1);
end

Gp = Gp.fit(connmatrixAll_D);
for i= 1:11
    pgPSI(i,:)=Gp.aligned{i}(:,1);
end

Gp = Gp.fit(connmatrixAll_M);
for i= 1:11
    pgMET(i,:)=Gp.aligned{i}(:,1);
end


plot_hemispheres([nanmean([-pg; -pg0;pg2]); nanmean(-pgMET); nanmean(pgPSI)]', ...
    {surf_lh,surf_rh},'parcellation',labeling, ...
    'labeltext',{'PLCB','MET','PSI'});

% Calculate gradient score difference

Gp = Gp.fit(connmatrixAll_B1);
for i = 1:length(connmatrixAll_B1)
    scoredif0(i) = max(Gp.aligned{i}(:,1))-min(Gp.aligned{i}(:,1));
end

Gp = Gp.fit(connmatrixAll_B2);
for i = 1:length(connmatrixAll_B2)
    scoredif(i) = max(Gp.aligned{i}(:,1))-min(Gp.aligned{i}(:,1));
end

Gp = Gp.fit(connmatrixAll_B3);
for i = 1:length(connmatrixAll_B3)
    scoredif2(i) = max(Gp.aligned{i}(:,1))-min(Gp.aligned{i}(:,1));
end

Gp = Gp.fit(connmatrixAll_D);
for i = 1:length(connmatrixAll_D)
    scoredifPSI(i) = max(Gp.aligned{i}(:,1))-min(Gp.aligned{i}(:,1));
end

Gp = Gp.fit(connmatrixAll_M);
for i = 1:length(connmatrixAll_M)
    scoredifMET(i) = max(Gp.aligned{i}(:,1))-min(Gp.aligned{i}(:,1));
end
