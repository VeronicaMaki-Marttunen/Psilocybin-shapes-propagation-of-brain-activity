% Estimates connectivity matrices

Subjects=[1:7 11 13 14 15];
%Subjects=[3,2,16,18,19,21,24,93,96,98,99];


tasks={'Baseline1','Baseline2','Baseline3','Drug','MET'};

runs=1;
num_bins=70;
Nsub=length(Subjects);
nsess=length(tasks);
nruns=1;

parcellation=ft_read_cifti('Schaefer2018_200Parcels_7Networks_order.dlabel.nii');

modules = parcellation.parcels(1:32492*2,:);
modulesn = length(unique(parcellation.parcels(~isnan(parcellation.parcels)))>0)-1;

connAll = nan(Nsub,nsess,modulesn,modulesn); % n subjs, n sessions


restingnt = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
for i = 1:length(restingnt)
  Index = find(contains(parcellation.parcelslabel,restingnt{i}));
modulesRSN(Index) = i;
end

for subject = 1:length(Subjects)

    disp(['Subject ',num2str(subject)])
    for sess = 1:nsess
disp(['Session ',num2str(sess)])

  inputfile = ['sub-',num2str(Subjects(subject)),'_',tasks{sess},'_rsfMRI_uout_bpss_sr_noGSR_sm4.dtseries.nii'];
       oi=ls (inputfile);
       
        if length(oi)>1
              [Mx_z_sparse,Mx_z_sparseGS] = calc_FC_function(inputfile,Subjects2(subject),modules,modulesn);
              connAll(subject,sess,:,:) = Mx_z_sparse;
              
        end
           
    end

end
save('conn_Schaef','connAll')

