function thresholded_matrix = calc_FC_function(inputfile,subj,modules,modulesn)

tr_data = 1.761; % TR
% Load regressors
if subj<10
    subjid=['0',num2str(subj)];
else
    subjid=num2str(subj);
end

cii_3 = ft_read_cifti(inputfile,'readdata',1);
dtser = cii_3.dtseries(1:32492*2,:);
nandtser = isnan(dtser(:,1));
dtser(nandtser,:) = [];

epi1msk = FT_Filter_mulch3(dtser',[0.001 .05],'bandpass',1/tr_data)'; % temporally smoothing data
epi1msk = zscore(epi1msk')';
gs_LR1 = mean(epi1msk); % calculate the global mean of input data

% Regress out GS - GSR
for i = 1:size(epi1msk,1)
    [B,BINT,R] = regress(squeeze(epi1msk(i,:))',squeeze(gs_LR1)');
    epi1msk_regGS(i,:) = R';
end


modules(nandtser,:) = [];

tmp_prin = epi1msk;

for mod1 = 1:modulesn
    temp = tmp_prin(modules == mod1,:);
    mod1tc = nanmean(temp);

    for mod2 = mod1+1:modulesn

        temp2 = tmp_prin(modules == mod2,:);
        mod2tc = nanmean(temp2);


        vx2 = [mod1tc ;mod2tc];
        a = corrcoef(vx2','Rows','complete');
        Mx(mod1,mod2)= a(2);

    end

end


% Z transformed
Mx_corr_z = atanh(Mx) ;
Mx_corr_z(end+1,:) = 0;
% Sparsity
numRows = size( Mx_corr_z, 1);
thresholded_matrix = zeros(numRows); % Initialize the output matrix

for i = 1:numRows
    % Get the z-scores for the current row
    row_scores =  Mx_corr_z(i, :);
    % Exclude self-connection
    row_scores(i) = -Inf; %
    % Determine the number of elements to retain
    num_to_keep = round(0.1 * numel(row_scores)); % 10% to keep (90% sparsity)
    % Get the indices of the strongest connections
    [~, sorted_indices] = sort(row_scores, 'descend');
    retain_indices = sorted_indices(1:num_to_keep);
    % Fill the thresholded matrix with the strongest connections
    thresholded_matrix(i, retain_indices) = row_scores(retain_indices);
end

clear mod1 mod2 Mx

% % For the results after regressing the global signal
% 
% tmp_prin = epi1msk_regGS;
% 
% for mod1 = 1:modulesn
%     temp = tmp_prin(modules == mod1,:);
%     mod1tc = nanmean(temp);
% 
%     for mod2 = mod1+1:modulesn
% 
%         temp2 = tmp_prin(modules == mod2,:);
%         mod2tc = nanmean(temp2);
% 
%         vx2 = [mod1tc ;mod2tc];
%         a = corrcoef(vx2','Rows','complete');
%         Mx(mod1,mod2)= a(2);
% 
%     end
% 
% end
% 
% 
% % Z transformed
% Mx_corr_z2 = atanh(Mx) ;
% Mx_corr_z2(end+1,:) = 0;
% % Sparsity
% numRows = size( Mx_corr_z2, 1);
% thresholded_matrix2 = zeros(numRows); % Initialize the 2nd output matrix
% 
% for i = 1:numRows
%     row_scores =  Mx_corr_z2(i, :);
%     row_scores(i) = -Inf; 
%     num_to_keep = round(0.1 * numel(row_scores)); % 10% to keep (90% sparsity)
%     [~, sorted_indices] = sort(row_scores, 'descend');
%     retain_indices = sorted_indices(1:num_to_keep);
%     thresholded_matrix2(i, retain_indices) = row_scores(retain_indices);
% end