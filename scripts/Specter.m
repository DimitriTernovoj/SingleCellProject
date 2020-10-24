% load all necessary library
addpath('Specter/dimred');
addpath('Specter/LSC');
addpath('Specter/utils');
addpath('Specter');
format short g;

%% read gene expression data: rows are cells, collumns are features (PCs, genes)
if downsampling == "0"
  data = csvread("../file_dir/specter_pca_data_normal.csv");
else
  data = csvread("../file_dir/specter_pca_data_sphetcher.csv");
end

specter_labels = eval_auto_Specter(data, n_clusters, ensemble_size, mingamma);

csvwrite("../file_dir/specter_clustering.csv", specter_labels);
