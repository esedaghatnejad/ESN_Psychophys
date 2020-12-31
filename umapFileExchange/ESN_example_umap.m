addpath('./umap/')

data_psort = Psort_extract_slot_data;
cs_wave = data_psort.topLevel_data.cs_wave;
window_ = 45:120;

[~, pca_mat, ~] = pca(cs_wave(:,window_));
[reduction, umap, clusterIdentifiers, extras]=run_umap(cs_wave(:,window_));

hFig = figure(1);
clf(hFig);

subplot(1, 3, 1)
plot(cs_wave')
xlabel('Time (ms)')
ylabel('Signal (uV)')

subplot(1, 3, 2)
plot(pca_mat(:, 1), pca_mat(:, 2), '.k')
xlabel('PCA1')
ylabel('PCA2')

subplot(1, 3, 3)
plot(reduction(:, 1), reduction(:, 2), '.k')
xlabel('UMAP1')
ylabel('UMAP2')

ESN_Beautify_Plot

