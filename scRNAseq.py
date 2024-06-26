import numpy as np
import pandas as pd
import os, time
import anndata
import scipy.stats
from itertools import compress
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib import rcParams
import re
import skmisc
import sys
import pickle
import scanpy as sc
import scanpy.external as sce
import copy
import cmasher as cmr
# import matplotlib.colors
import matplotlib.colors as clr
import sys
import session_info
sys.setrecursionlimit(100000)

session_info.show()
EPS = 1


# all this class used the package Scanpy- single-cell RNA sequencing analysis
class scRNAseqTA:
    def __init__(self, info_file,  all_process=True):
        # figures scanpy setting
        sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
        sc.logging.print_versions()
        sc.settings.set_figure_params(dpi=80, facecolor='white') #   , vector_friendly=False )
        sc.settings.autosave = True
        sc.settings.autoshow = False
        self.date = time.strftime('%Y%m%d_%H%M')  # time.strftime('_%Y-%m-%d_%H:%M')
        self.convert = pd.read_csv('genes_names7.csv')
        self.mutants = pd.read_csv('FishStrains.csv') #  markersTest
        self.info = pd.read_csv(info_file.replace('txt', 'csv'))
        self.org = 'kf'
        self.DEgenes = []
        # print('make sure you filled correct the information experiment file.\nYou can fill manual or use: fill_scRNAseq_info.py')
        self.read_info_file(info_file)
        self.read_files_concatenate()
        # self.diffusion_pseudotime()
        # self.filtering_normalization_scaling()
        # self.pca_and_umap(res=0.3)
        # self.cell_type_identification('leiden')

    def read_info_file(self, info_file):
        """
        The function takes an info file with specific pattern and save the parameters
        the parameters include: path to data folder, names of folderName, group info
        :param info_file: an info file with detail of the folder and saved data
        :return:
        """
        with open(info_file, 'r') as inf:
            field = inf.read().split('\n')
            self.path = field[2].split(':')[1].strip()
            self.folderName = field[3].split(':')[1].strip().replace(' ', '').split(',')
            self.analysis_path = field[4].split(':')[1].strip()
            os.popen('mkdir -p ' + self.analysis_path + self.date)
            sc.settings.figdir = self.analysis_path + self.date

            parameters = field[6].split(';')
            up_to_word = ":"
            rx_to_first = r'^.*?{}'.format(re.escape(up_to_word))
            self.parameters = dict(zip(['count_min', 'count_max', 'highVG', 'neighbors', 'pcs', 'resolution'],
                                       [float(re.sub(rx_to_first, '', x, flags=re.DOTALL).strip()) for x in parameters]))
            os.popen('cp ' + info_file + ' ' + self.analysis_path + self.date)

            print(self.parameters)

    def read_files_concatenate(self):
        """
        The function goes over the expression files of each Group table (cells in rows, genes in columns) and concatenate them.
        Then, merge the info table (with batch and group information for each experiment) in the observation table (.obs)
        Finally, save sc matrix in self parameters
        :return:
        """
        ExpMatrix = {}
        i = 0
        for b in self.folderName:
            matrix_name = 'named_RSEM_%s_matrix_expected_counts_%s.csv' % (self.org, b)
            # matrix_name = 'RSEM_NFU_matrix_UNIQUE_counts_%s.bed' % b
            exM = sc.read_csv(os.path.join(self.path, matrix_name), delimiter='\t', first_column_names=True)
            exM.var_names_make_unique()
            exM.obs_names_make_unique()
            sc.pp.filter_cells(exM, min_genes=1)
            sc.pp.filter_genes(exM, min_cells=1)


            exM.obs_names = exM.obs_names + str(i)
            i += 1
            print(exM)
            ExpMatrix[b] = exM


        self.sc_matrix = anndata.concat(ExpMatrix, label="mergedGroup", join='outer', fill_value=0)
        print(self.sc_matrix.obs)
        self.sc_matrix.obs = pd.merge(self.sc_matrix.obs, self.info, left_on='mergedGroup', right_on='FolderName', how='left').set_axis(self.sc_matrix.obs.index)
        # self.sc_matrix.obs['Batch'] = self.sc_matrix.obs['Batch'].astype('category')
        print(self.sc_matrix.obs)
        # sc.pp.combat(self.sc_matrix)

        # if self.num_groups != len(self.groups):
        #     print(self.groups)
        # self.sc_matrix.obs = self.sc_matrix.obs.replace(self.groups)
        self.sc_matrix.var_names = [g.replace('gene-', '') for g in self.sc_matrix.var_names]
        self.sc_matrix.var_names = self.__convertGenesNames(self.sc_matrix.var_names) # converting the gene name to NCBI final symbol

    def calculate_properties_on_group(self, g1, g2='rest', cluster='Group'):
        """
        The function adding the average expression of the cells in the indicated group according self.group to genes annotation(.var)
        In addition, the function calculate log2 fold change between g1 qnd g2
        The function use the raw values

        :param g1: group name
        :param g2: group name, by deflate: 'rest'- calculate the average expression of all the rest cells
        :param cluster:
        :return:
        """
        self.sc_matrix.var[g1 + '_observation'] = np.mean(np.expm1(self.sc_matrix.raw[self.sc_matrix.obs[cluster] == g1, :].X), axis=0)
        if g2 == 'rest':
            self.sc_matrix.var['rest_observation'] = np.mean(np.expm1(self.sc_matrix.raw[self.sc_matrix.obs[cluster] != g1, :].X), axis=0)
        else:
            self.sc_matrix.var[g2 + '_observation'] = np.mean(np.expm1(self.sc_matrix.raw[self.sc_matrix.obs[cluster] == g2, :].X), axis=0)
        print(self.sc_matrix.var)

        self.sc_matrix.var['total'] = np.count_nonzero(self.sc_matrix.X, axis=0)
        self.sc_matrix.var['logFC'] = np.log2(self.sc_matrix.var[g1 + '_observation']+EPS) -  np.log2(self.sc_matrix.var[g2 + '_observation']+ EPS)

    def __convertGenesNames(self, geneList, source='NCBI', target='Final symbol'):
        """
        The function convert gene names between 1.human to symbol names, 2.NCBI names, 3.final symbol of the killifish
        source and target could be the column name of the given table (self.convert). in this case: NCBI / Final symbol/ Human
        the source should be identical to the gene annotation (deflate- NCBI)
        :param geneList: list of the old name
        :param source: column name of the source gene name
        :param target: column name of the source gene name
        :return: list of the converted gene names
        """
        if type(geneList) is not list:
            geneList = geneList.to_list()
        return self.convert.set_index(source).loc[geneList].reset_index()[target].to_list()

    def filtering_normalization_scaling(self):
        """
        filtering un-informative cells and genes, and throw out cell with extreme expression
        normalization the data
        scaling to normal distribution

        The parameters are from the info file
        it is recommend first use default values and afterward chose the correct values from violin plot
        :return:
        """
        sc.pl.highest_expr_genes(self.sc_matrix, n_top=20, save=self.date, show=False, )
        sc.pp.filter_cells(self.sc_matrix, min_genes=1)
        sc.pp.filter_genes(self.sc_matrix, min_cells=1)
        self.sc_matrix.obs['n_counts'] = self.sc_matrix.X.sum(axis=1)

        self.sc_matrix.var['mt'] = self.sc_matrix.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(self.sc_matrix, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        print('*****************************************')
        print(self.sc_matrix)

        sc.pl.violin(self.sc_matrix, ['n_genes', 'n_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save='before_'+self.date, show=False)
        sc.pl.violin(self.sc_matrix, ['n_genes', 'n_counts', 'pct_counts_mt'], groupby='FolderName', jitter=0.4, multi_panel=True, save='group_'+self.date, rotation=90, show=False)
        print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
        print(pd.crosstab(index=self.sc_matrix.obs['Group'], columns='count'))

        sc.pl.scatter(self.sc_matrix, x='n_genes', y='n_counts', save='_before', show=False)

        self.sc_matrix = self.sc_matrix[self.parameters['count_max'] > self.sc_matrix.obs.n_genes, :]
        self.sc_matrix = self.sc_matrix[self.parameters['count_min'] < self.sc_matrix.obs.n_genes, :]
        self.sc_matrix = self.sc_matrix[self.sc_matrix.obs.pct_counts_mt < 20, :]

        sc.pl.violin(self.sc_matrix, ['n_genes', 'n_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save=self.date, show=False)
        sc.pl.violin(self.sc_matrix, ['n_genes', 'n_counts', 'pct_counts_mt'], groupby='FolderName', jitter=0.4, multi_panel=True, save='_groupAfter_'+self.date, rotation=90, show=False)
        sc.pl.scatter(self.sc_matrix, x='n_genes', y='n_counts', save='_after', show=False)

        print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
        print(pd.crosstab(index=self.sc_matrix.obs['Group'], columns='count'))

        self.sc_matrix.obs['200'] = (np.count_nonzero(self.sc_matrix.X, axis=1) < 200).astype(int)
        print(self.sc_matrix.obs)
        print(np.mean(self.sc_matrix.obs.n_counts), np.median(self.sc_matrix.obs.n_counts))
        print(np.mean(self.sc_matrix.obs.n_genes), np.median(self.sc_matrix.obs.n_genes))

        sc.pp.highly_variable_genes(self.sc_matrix, flavor='seurat_v3', n_top_genes=int(self.parameters['highVG'])) #) #min_mean=0.0125, max_mean=3, min_disp=0.5,
        HV = self.sc_matrix.var
        HV['NCBI'] = self.__convertGenesNames(HV.index, source='Final symbol', target='NCBI')
        self.mutants = self.mutants[self.mutants['FinalSymbol'].isin(self.sc_matrix.var_names)]['FinalSymbol'].tolist()

        # self.diffusion_pseudotime()
        print(self.sc_matrix)

        # normalize
        sc.pp.normalize_total(self.sc_matrix, target_sum=1e4)
        print("%%%%%%%%%", sum(sum(self.sc_matrix.X)))
        sc.pp.log1p(self.sc_matrix)
        print("%%%%%%%%%", sum(sum(self.sc_matrix.X)))

        self.sc_matrix.raw = self.sc_matrix

        # scaling
        sc.pp.scale(self.sc_matrix, max_value=10)
        print("%%%%%%%%%", sum(sum(abs(self.sc_matrix.X))))
        print('#############################################')
        print(self.sc_matrix)
        print('==========================================================')
        print(pd.crosstab(index=self.sc_matrix.obs['Group'], columns='count'))

    def pick_specific_gene_list(self, gene_list=[]):
        """
        The function subset the given gene list

        :param gene_list: list of genes that required in this analysis
        if gene list is empty - the function returned the whole sc matrix
        :return: sc matrix with gene list
        """
        if not len(gene_list):
            sc_matrix = self.sc_matrix
        else:
            sc_matrix = self.sc_matrix[:, self.sc_matrix.var_names.isin(gene_list)]
        return sc_matrix

    def pick_specific_cells(self, cluster='Group', clusterNamber=[]):
        """
        The function subset cells in the indicated cluster

        :param cluster: the name of the cluster to filter accordingly
        :param clusterNamber: the number of the cluster
        :return:
        """
        if len(clusterNamber) == 0:
            return self.sc_matrix
        else:
            print(self.sc_matrix.obs[cluster].isin(clusterNamber))
            sc_matrix = self.sc_matrix[self.sc_matrix.obs[cluster].isin(clusterNamber), :]
        return sc_matrix

    def pick_specific_cells_by_name(self, cellList):
        """
        The function subset the given cell list
        :param cellList: list of the cell
        :return:
        """
        sc_matrix = self.sc_matrix[self.sc_matrix.obs_names.isin(cellList), :]
        return sc_matrix

    def pca_and_umap(self, gene_list=[], neighbors=0.0, pcs=0.0, res=0.0, dist=0.5, repSource=None, repTarget=None, name='', file_name='', *args, **kwargs):
        """
        The function calculate and draw PCA and UMAP.

        :param gene_list: list of interesting genes
        :param neighbors: number of neighbors
        :param pcs: number of PCs
        :param res: cluster resolution
        :param dist: The effective minimum distance between embedded points. (sc.tl.umap())

        if you want to change the cluster name put lists in the next two parameters. they should be with same length
        :param repSource: list of cluster names (could be list of list, while you group several cluster under one name)
        :param repTarget: list of target cluster name

        :param name: add the string to the name of each plot
        :param args: parameter for (sc.pl.umap())
        :param kwargs: parameter for (sc.pl.umap())
        :return: --
        """

        # using the parameter from the info file by deflate
        if neighbors == 0.0:
            neighbors = self.parameters['neighbors']
        if pcs == 0.0:
            pcs = self.parameters['pcs']
        if res == 0.0:
            res = self.parameters['resolution']

        sc_matrix = self.pick_specific_gene_list(gene_list)
        # PCA
        sc.tl.pca(sc_matrix, svd_solver='arpack')
        sc.pl.pca(sc_matrix, color='Group', components=['1,2', '3,4', '5,6', '7,8'], ncols=2, save=self.date, show=False)
        sc.pl.pca_variance_ratio(sc_matrix, log=True, save=self.date , show=False)
        sce.pp.harmony_integrate(sc_matrix, 'Batch')  # batch correction
        sc.pp.neighbors(sc_matrix, n_neighbors=int(neighbors), n_pcs=int(pcs), use_rep='X_pca_harmony') #)

        #
        sc.pp.neighbors(sc_matrix, n_neighbors=int(neighbors), n_pcs=int(pcs), use_rep='X_pca_harmony') #
        # UMAP
        sc.tl.leiden(sc_matrix, resolution=res) #
        sc.tl.louvain(sc_matrix, resolution=res)

        star = ['star (1 of 2)', 'star (2 of 2)', 'cyp17a1 (1 of 2)', 'cyp17a1 (2 of 2)', 'cyp11a1 (1 of 2)', 'cyp11a1 (2 of 2)', 'ddx4']
        sc.tl.umap(sc_matrix, min_dist=dist)

        self.sc_matrix.obs['cellType'] = self.sc_matrix.obs['leiden']
        print(pd.crosstab(index=self.sc_matrix.obs['leiden'], columns='count'))

        if (repSource is not None) & (repTarget is not None):
            self.sc_matrix.obs[['cellType']].replace(repSource, repTarget, inplace=True)

        sc.pl.umap(sc_matrix, color=['leiden', 'louvain', 'Group', 'Batch', 'n_genes', 'n_counts', 'mergedGroup', 'cellType'], show=False, save=self.date + name)
        sc.pl.umap(sc_matrix, color=['Group', 'cellType'], show=False,  save=self.date + 'size' + name, *args, **kwargs)

        for g in list(set(self.sc_matrix.obs['Group'])):
            sc.pl.umap(sc_matrix, color=['Group'], groups=g, show=False, save=g+name)
        for b in list(set(self.sc_matrix.obs['Batch'])):
            sc.pl.umap(sc_matrix, color=['Batch'], groups=b, show=False, save=str(b) + name)

        sc.pl.umap(sc_matrix, color=self.mutants, use_raw=False, show=False, save='_OurFishStrains' + name)

        if len(file_name) > 0:
            genes = pd.read_csv(file_name)['Final Symbol'].to_list()
            genes = [m for m in genes if m in self.sc_matrix.var_names]

            sc.pl.umap(sc_matrix, color=genes, use_raw=False, show=False, save='_specific_markers' + name, *args, **kwargs)

    def differential_expression_genes(self, g1, g2='rest', method='wilcoxon', cluster='Group', name='', thresholdpval=0.05, thresholdFC=1.0):
        """
        The function calculate the differential gene expression between g1 and g2.
        g1 and g2 are cluster number from the parameter cluster

        :param g1:
        :param g2:
        :param method: statistics test
        :param cluster:
        :param name: add the string to the name of each plot
        :param thresholdpval: threshold of the p-value
        :param thresholdFC: threshold of the FC
        :return:
        """
        params = ['names', 'scores', 'pvals', 'pvals_adj']#, 'pts', 'pts_rest']

        # to filter the genes with no expression in further analysis
        if g2 == 'rest':
            self.sc_matrix.var['maxEx'] = np.max(self.sc_matrix.raw.X, axis=0)
            self.sc_matrix.var['meanEx'] = np.mean(self.sc_matrix.raw.X, axis=0)
        else:
            self.sc_matrix.var['maxEx'] = np.max(self.sc_matrix.raw[self.sc_matrix.obs[cluster].isin([g1, g2]), :].X, axis=0)
            self.sc_matrix.var['meanEx'] = np.mean(self.sc_matrix.raw[self.sc_matrix.obs[cluster].isin([g1, g2]), :].X, axis=0)

        # compare two clusters
        sc.tl.rank_genes_groups(self.sc_matrix, cluster, groups=[g1], reference=g2, corr_method='bonferroni', method=method, pts=True, tie_correct=True) #, n_genes=-1)
        print(self.sc_matrix.uns['rank_genes_groups'])
        summary_table = pd.DataFrame([self.sc_matrix.uns['rank_genes_groups'][x][g1] for x in params], index=params).T
        self.calculate_properties_on_group(g1, g2, cluster=cluster) # calculate the log2 FC (in the way Seurat do)
        significant_genes_table = summary_table.merge(self.sc_matrix.var[['logFC', 'maxEx', 'meanEx']], left_on='names', right_index=True)
        significant_genes_table['NCBI'] = self.__convertGenesNames(significant_genes_table['names'],source='Final symbol', target='NCBI')
        significant_genes_table['Human'] = self.__convertGenesNames(significant_genes_table['names'],source='Final symbol', target='Human')
        significant_genes_table.sort_values('logFC', inplace=True)
        significant_genes_table.to_csv(os.path.join(self.analysis_path, self.date, '%s.csv' % name))

        genes = significant_genes_table[(abs(significant_genes_table['logFC']) > thresholdFC) & (significant_genes_table['pvals_adj'] < thresholdpval)]
        genes.sort_values(by=['logFC'], inplace=True)
        self.DEgenes.extend(genes['names'].to_list())

        sigGenes = significant_genes_table.head(10)['names'].to_list()
        sigGenes.extend(significant_genes_table.tail(10)['names'].to_list())
        sc.pl.stacked_violin(self.sc_matrix, sigGenes, groupby=cluster, rotation=90, show=False, save=name)

    def marker_genes(self, gene_list=[], method='wilcoxon', cluster='Group', name='', thresholdpval=0.05, thresholdFC=1.0):
        """
        The function calculate the DGE between each cluster to the remaining cluster

        :param gene_list: specific gene list to look at
        :param method: statistic test. options: ['t-test', '', 't-test_overestim_var','logreg'] default: 'wilcoxon'
        :param name: add the string to the name of each plot
        :param thresholdpval: threshold of the p-value
        :param thresholdFC: threshold of the FC
        :return: dictionary of significant gene list in each group. keys-group, value-gene list
        """
        params = ['names', 'scores', 'pvals', 'pvals_adj', 'pts', 'pts_rest']#]
        sc_matrix = self.pick_specific_gene_list(gene_list)

        # compare each cluster to the rest
        marker_genes = {}

        sc.tl.rank_genes_groups(sc_matrix, cluster, corr_method='bonferroni', method=method, pts=True, tie_correct=True) #, n_genes=-1)
        # sc.tl.rank_genes_groups(sc_matrix, cluster, groups=[0], reference='1' , method=method, n_genes=-1)
        print('**************************************************')
        sc.pl.rank_genes_groups(sc_matrix, show=False, save=self.date)
        print(sc_matrix.uns['rank_genes_groups'])
        for g in set(self.sc_matrix.obs[cluster]):
            # calculate the p-value and score by t-test
            summary_table = pd.DataFrame([sc_matrix.uns['rank_genes_groups'][x][g] for x in params], index=params).T
            self.calculate_properties_on_group(g, cluster=cluster)
            significant_genes_table = summary_table.merge(self.sc_matrix.var[['logFC']], left_on= 'names', right_index=True)
            significant_genes_table['NCBI'] = self.__convertGenesNames(significant_genes_table['names'], source='Final symbol', target='NCBI')
            significant_genes_table['Human'] = self.__convertGenesNames(significant_genes_table['names'],source='Final symbol', target='Human')
            significant_genes_table['NCBI Definition'] = self.__convertGenesNames(significant_genes_table['names'],source='Final symbol', target='NCBI Definition')
            marker_genes[g] = significant_genes_table
            significant_genes_table.to_csv(os.path.join(self.analysis_path, self.date,  name + g + cluster + '.csv'))

            genes = significant_genes_table[(abs(significant_genes_table['logFC']) > thresholdFC) & (significant_genes_table['pvals_adj'] < thresholdpval)]
            genes.sort_values(by=['logFC'], inplace=True)
            genes.to_csv(os.path.join(self.analysis_path, self.date, 'significant_differential_gene_exp_up_%s%s_pval_%.2f_FC_%.1f.csv' % (name, g, thresholdpval, thresholdFC)))
           
        sc.pl.rank_genes_groups_heatmap(self.sc_matrix, n_genes=20, groupby=cluster, show_gene_labels=True, save=name + '1')
        sc.pl.rank_genes_groups_heatmap(self.sc_matrix, n_genes=20, groupby='Group', save=name + '2')
        return marker_genes

    def diffusion_pseudotime(self, name=''):
        """
        note: this function works on original data, without normalization or scaling

        :return:
        """
        
        # sc.pp.recipe_zheng17(self.sc_matrix)
        sc.pp.neighbors(self.sc_matrix, n_neighbors=20, method='gauss', use_rep='X_pca_harmony')
        sc.tl.diffmap(self.sc_matrix, n_comps=10)
        # sc.pp.neighbors(self.sc_matrix, n_neighbors=20, use_rep='X_diffmap')
        self.sc_matrix.obs['irootCluster'] = 'A'
        print(self.sc_matrix.obs['irootCluster'])
        self.sc_matrix.obs['irootCluster'][np.flatnonzero(self.sc_matrix.obs['leiden'] == '2')[117]] = 'B'
        print(self.sc_matrix.obs['irootCluster'][np.array(self.sc_matrix[:, 'piwil1'].X).argmax()])
        self.sc_matrix.uns['iroot'] = np.flatnonzero(self.sc_matrix.obs['leiden'] == '2')[117]      #np.flatnonzero(self.sc_matrix.obs['leiden'] == '2')[5] np.array(self.sc_matrix[:, 'piwil1'].X).argmax()

        sc.tl.dpt(self.sc_matrix, n_branchings=1)
        sc.pl.dpt_groups_pseudotime(self.sc_matrix, show=False, save='Pseudotime' + name + self.date)
        sc.pl.diffmap(self.sc_matrix, components=['1,2', '1,3', '2,3', '2,4'], color=['dpt_pseudotime', 'dpt_groups','leiden', 'Group'], show=False, save='Group' + name)

        sc.tl.umap(self.sc_matrix)
        sc.pl.umap(self.sc_matrix, color=['Group', 'cellType', 'dpt_pseudotime', 'n_genes', 'n_counts', 'mergedGroup', 'irootCluster'], show=False,  save=self.date + name + 'diff')

        aa = self.sc_matrix.obs.groupby(['leiden']).mean()  # 'dpt_pseudotime'], )
        print(aa['dpt_pseudotime'])
        #
        print(np.array(self.sc_matrix[:, 'piwil1'].X))
        print(np.array(self.sc_matrix[:, 'piwil1'].X).argmax())

    def cell_type_identification(self, clusterMethod='leiden', name='',markerGeneFile='all', order=[], *args, **kwargs):
        """
        The function draw dot-plot of the provided gene-list

        :param clusterMethod: name of the cluster method
        :param name: add the string to the name of each plot
        :param markerGeneFile: file name for the plot, if deflate- print all the options
        :param order: the order of the cluster names
        :param args:
        :param kwargs:
        :return:
        """
        cmap = cmr.get_sub_cmap('PuBu', 0.1, 1.0)
        # cmap = clr.LinearSegmentedColormap.from_list('3', ['#E6E6E6', '#039DDB', '#006497', '#022234'], N=256)
        cmap = clr.LinearSegmentedColormap.from_list('3', ['#E6E6E6', '#54B7DF', '#006497', '#022234'], N=256)

        if len(order) > 0:
            # print(scscWT.sc_matrix.obs['GermClusters'])
            self.sc_matrix.obs[clusterMethod].cat.set_categories(order, inplace=True)
            # print(scscWT.sc_matrix.obs['GermClusters'])

        sc.tl.dendrogram(self.sc_matrix, groupby=clusterMethod)

        sc.tl.rank_genes_groups(self.sc_matrix, groupby=clusterMethod, corr_method='bonferroni', method='wilcoxon', pts=True, tie_correct=True)
        sc.pl.rank_genes_groups_dotplot(self.sc_matrix, n_genes=10, save=self.date + name)  # , groups=groupDots)

        if markerGeneFile =='all':
            # Testis PMID: 30315278
            markerGenes = pd.read_csv(
                '/cellType/markersHumanPaperTestis2018.csv')
            markerHpaper = markerGenes.groupby('cluster')['Human'].apply(list).to_dict()
            for h in markerHpaper:
                markerHpaper[h] = [m for m in markerHpaper[h] if m in self.convert['Human'].tolist()]
                markerHpaper[h] = [m for m in self.__convertGenesNames(markerHpaper[h], source='Human') if
                                   m in self.sc_matrix.var_names]
            print(markerHpaper)
            sc.pl.dotplot(self.sc_matrix, markerHpaper, clusterMethod, dendrogram=True, save='_Testis' + self.date + name)
            # human papers
            markerGenes = pd.read_csv(
                '/cellType/markersHumanPaper.csv')  #Ovarian2020
            markerHpaper = markerGenes.groupby('cluster')['Human'].apply(list).to_dict()
            for h in markerHpaper:
                markerHpaper[h] = [m for m in markerHpaper[h] if m in self.convert['Human'].tolist()]
                markerHpaper[h] = [m for m in self.__convertGenesNames(markerHpaper[h], source='Human') if
                                   m in self.sc_matrix.var_names]
            print(markerHpaper)
            sc.pl.dotplot(self.sc_matrix, markerHpaper, clusterMethod, dendrogram=True, save='_H_paper' + self.date + name)

            markerGenes = pd.read_csv('/cellType/markersHumanTestis.csv') #markersHumanPaperTestis2018
            markerHpaper = markerGenes.groupby('cluster')['Human'].apply(list).to_dict()
            for h in markerHpaper:
                markerHpaper[h] = [m for m in markerHpaper[h] if m in self.convert['Human'].tolist()]
                markerHpaper[h] = [m for m in self.__convertGenesNames(markerHpaper[h], source='Human') if m in self.sc_matrix.var_names]
            print(markerHpaper)
            sc.pl.dotplot(self.sc_matrix, markerHpaper, clusterMethod, dendrogram=True, save='_H_testis' + self.date + name)

            # Human protein atlas
            markerGenes = pd.read_csv(
                '/cellType/markersHumanProteinAtlas.csv')
            markerHpaper = markerGenes.groupby('cluster')['Human'].apply(list).to_dict()
            for h in markerHpaper:
                markerHpaper[h] = [m for m in markerHpaper[h] if m in self.convert['Human'].tolist()]
                markerHpaper[h] = [m for m in self.__convertGenesNames(markerHpaper[h], source='Human') if
                                   m in self.sc_matrix.var_names]
            print(markerHpaper)
            sc.pl.dotplot(self.sc_matrix, markerHpaper, clusterMethod, dendrogram=True, save='_HumanAtlas' + self.date + name)

            # immune killifish
            markerGenes = pd.read_csv(
                '/cellType/immune_markers.csv')
            markerHpaper = markerGenes.groupby('cluster')['NCBI'].apply(list).to_dict()
            for h in markerHpaper:
                markerHpaper[h] = [m for m in markerHpaper[h] if m in self.convert['NCBI'].tolist()]
                markerHpaper[h] = [m for m in self.__convertGenesNames(markerHpaper[h], source='NCBI') if
                                   m in self.sc_matrix.var_names]
            print(markerHpaper)
            sc.pl.dotplot(self.sc_matrix, markerHpaper, clusterMethod, dendrogram=True, save='_Immune' + self.date + name)

            # zebrafish
            # the DE genes for each cell type by PMID: 35588359
            with open('/cellType/dictCellTypesZFgenes.pickle',
                      'rb') as handle:
                markerZFpaper = pickle.load(handle)

            for m in markerZFpaper:
                markerZFpaper[m] = [m for m in markerZFpaper[m] if m in self.sc_matrix.var_names]
            sc.pl.dotplot(self.sc_matrix, markerZFpaper, clusterMethod, dendrogram=True, save='_ZF_' + self.date + name)

        else:
            # in case you provide your own gene list
            markerGenes = pd.read_csv('%s.csv' % markerGeneFile)
            markerK = markerGenes['Final Symbol'].to_list()
            markerK = [m for m in markerK if m in self.sc_matrix.var_names]
            print(markerK)
            sc.pl.dotplot(self.sc_matrix, markerK, groupby=clusterMethod, cmap='PuBu', save='_MyK_PuBu' + self.date + name)
            sc.pl.dotplot(self.sc_matrix, markerK, groupby=clusterMethod, cmap=cmap, save='_MyK_cmap' + self.date + name)
            sc.pl.umap(self.sc_matrix, color=['Group', clusterMethod], show=False, save= '_' + clusterMethod + name, *args, **kwargs)
            sc.pl.umap(self.sc_matrix, color=markerK, show=False, color_map = cmap, save='_MyK'+ name, *args, **kwargs)
            sc.pl.umap(self.sc_matrix, color=markerK, show=False, save='_MyK2'+ name, *args, **kwargs)

        # All markers
        markerGenes = pd.read_csv('markersL.csv')
        markerK = markerGenes.groupby('cluster')['Final Symbol'].apply(list).to_dict()
        for k in markerK:
            markerK[k] = [m for m in markerK[k] if m in self.sc_matrix.var_names]

        sc.pl.dotplot(self.sc_matrix, markerK, clusterMethod, dendrogram=True, save='_MyL' + self.date + name)
        genes = [m for m in markerGenes['Final Symbol'] if m in self.sc_matrix.var_names]
        sc.pl.umap(self.sc_matrix, color=genes, show=False, save='_MyL' + name, *args, **kwargs)

        print("===================================================")
        print(pd.crosstab(index=self.sc_matrix.obs[clusterMethod], columns='count'))
        print(pd.crosstab(index=self.sc_matrix.obs['mergedGroup'], columns='count'))
        print(pd.crosstab(index=self.sc_matrix.obs['Batch'], columns='count'))

    def FCheatmap(self, geneList, cmap='bwr', name='matrix'):
        """
        The function draw heatmap of the FC between two following columns.
        use the function after "differential_expression_genes"

        :param geneList: gene list of interesting
        :param cmap:
        :param name: add the string to the name of each plot
        :return: --
        """
        matrix = self.sc_matrix.var.loc[geneList, :]
        matrix = matrix[[col for col in matrix.columns if '_observation' in col]]
        matrix.columns = [t.replace('_observation', '') for t in matrix.columns.to_list()]
        fc = matrix.copy()
        for i in range(0, fc.shape[1], 2):
            fc.iloc[:, i:i + 2] = fc.iloc[:, i:i + 2].div(fc.iloc[:, i].to_list(), 0)
        fc.replace([np.inf, -np.inf], np.nan, inplace=True)
        fc.replace([np.nan], 1, inplace=True)
        LOGfc = np.log2(fc).iloc[:,[i for i in range(1,fc.shape[1], 2)]]
        print(LOGfc)
        divnorm = mpl.colors.TwoSlopeNorm(vmin=-abs(LOGfc).max().max(), vcenter=0., vmax=abs(LOGfc).max().max())
        fig, ax = plt.subplots()
        ax.imshow(LOGfc, cmap=cmap, norm=divnorm)
        ax.grid(False)
        ax.set_yticks(np.arange(LOGfc.shape[0]), labels=LOGfc.index)
        ax.set_xticks(np.arange(LOGfc.shape[1]), labels=LOGfc.columns)
        plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")
        fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=divnorm))
        # plt.show()
        plt.savefig(self.analysis_path + self.date + '/' + name + '.pdf')

    def subsampling_by_group(self, groupSub='KOF', groupFilterby='WTF'):
        scscCopy = copy.deepcopy(self)
        scscCopy.sc_matrix = self.pick_specific_cells('Group', [groupSub])
        sc.pp.subsample(scscCopy.sc_matrix, n_obs=int(pd.crosstab(index=self.sc_matrix.obs['Group'], columns='count').loc[groupFilterby]))  #
        print(len(list(scscCopy.sc_matrix.obs_names)))
        print(len(list(self.sc_matrix[self.sc_matrix.obs['Group'] == groupFilterby, :].obs_names)))
        return self.pick_specific_cells_by_name(list(scscCopy.sc_matrix.obs_names) + list(self.sc_matrix[self.sc_matrix.obs['Group'] == groupFilterby, :].obs_names))

scsc = scRNAseqTA("info.txt")


# ##################WT################
scscWT = copy.deepcopy(scsc)
scscWT.sc_matrix = scsc.pick_specific_cells('Group', ['WTM', 'WTF'])
scscWT.filtering_normalization_scaling()
a = [['0', '1', '2', '3', '4', '12'], ['13'], ['5', '6', '7', '8', '9', '10', '11', '14', '15']]
b=['Germ', 'Eggs', 'Somatic']

scscWT.pca_and_umap(neighbors=40.0, pcs=20.0, res=0.6, dist=0.5, size=15, repSource=a, repTarget=b)
sc.pl.violin(scscWT.sc_matrix, ['n_genes', 'n_counts', 'pct_counts_mt'], groupby='cellType', jitter=0.4, multi_panel=True, save='_cellTypeAfter_', rotation=90, show=False)

scscWT.cell_type_identification()
scscWT.cell_type_identification(size=15, markerGeneFile='markersGerm', name='GermALL')

scscWT.sc_matrix.obs['leiden2'] = scscWT.sc_matrix.obs['leiden']
scscWT.sc_matrix.obs[['GermClusters']].replace(a[2], ['Somatic'], inplace=True)

scscWT.sc_matrix.obs['GermClusters'] = (scscWT.sc_matrix.obs['leiden2'].astype(str) + scscWT.sc_matrix.obs['Group'].astype(str)).astype('category')
d = [['3WTM', '5WTM', '6WTM', '8WTM'], ['3WTF', '5WTF', '6WTF', '7WTF', '10WTF'], ['2WTF', '1WTF']]
scscWT.sc_matrix.obs[['GermClusters']].replace(d, ['Somatic', 'Somatic', '1WTF'], inplace=True)


newOrder = ['9WTF', '1WTF', '0WTF', '2WTM', '4WTM', '1WTM', '0WTM', 'Somatic']

# ########### WT somatic ###############
scscWTsom = copy.deepcopy(scscWT)
scscWTsom.sc_matrix = scscWT.pick_specific_cells('cellType', ['Somatic'])
scscWTsom.pca_and_umap(neighbors=20.0, pcs=15.0, res=0.4, name='WTSomatic', size=30)


# scscWTsom.sc_matrix = scscWTsom.subsampling_by_group('WTF', 'WTM')

print(pd.crosstab(index=scscWTsom.sc_matrix.obs['Group'], columns='count'))

# newOrder = ['11', '9', '5', '7', '2', '8', '4', '1', '6', '0', '10', '3'], order=newOrder
newOrder = ['3', '11', '0', '7', '1', '4', '9', '2', '6', '5', '10', '8']
scscWTsom.cell_type_identification(name='WTSomatic', size=30, order=newOrder, markerGeneFile='markersWTcellTypes')  #)
newOrder = ['8', '7', '0', '10', '5', '6', '2', '4', '1', '9', '11', '3']
# scscWTsom.cell_type_identification(name='WTSomaticZF', size=30, order=newOrder, markerGeneFile='markersCompareZF')
# scscWTsom.marker_genes(cluster='leiden', name='WTSomatic')


# ########### WT Germ ###############
# WT Germ (male)
scscWTg = copy.deepcopy(scscWT)
scscWTg.sc_matrix = scscWT.pick_specific_cells('cellType', ['Germ', 'Eggs'])
scscWTg.pca_and_umap(name='Germ', neighbors=40.0, pcs=20.0, res=0.6, dist=0.4)   #,0.6
scscWTg.sc_matrix.obs['GermClusters'] = (scscWTg.sc_matrix.obs['leiden'].astype(str) + scscWTg.sc_matrix.obs['Group'].astype(str)).astype('category')  #
# scscWTg.diffusion_pseudotime()
# scscWTg.cell_type_identification(name='Germ', markerGeneFile='markersGerm')  #)

d = [['3WTF', '0WTF'], ['6WTF', '1WTF']]
scscWTg.sc_matrix.obs[['GermClusters']].replace(d, ['0WTF', '1WTF'], inplace=True)
sc.pl.umap(scscWTg.sc_matrix, color=['Group', 'GermClusters'], show=False, save=scscWTg.date + 'sizeGerm2')

print(pd.crosstab(index=scscWTg.sc_matrix.obs['Group'], columns='count'))
newOrder = ['6WTM',     '8WTF','1WTF','0WTF','2WTF',        '3WTM','4WTM','7WTM','1WTM','0WTM','5WTM','2WTM']
scscWTg.cell_type_identification(clusterMethod='GermClusters', name='GermClusters', order=newOrder, markerGeneFile='markersGerm') #)

scscWTg2 = copy.deepcopy(scscWTg)
scscWTg2.sc_matrix = scscWTg.pick_specific_cells('leiden', ['0', '1', '2', '3', '4', '5', '7'])
scscWTg2.sc_matrix = scscWTg2.pick_specific_cells('Group', ['WTM'])
scscWTg2.pca_and_umap(name='Germ2', neighbors=40.0, pcs=20.0, res=0.6, dist=0.4)   #,0.6

# scscWTg2.marker_genes(cluster='leiden', name='WTgerm')
scscWTg2.diffusion_pseudotime(name='Germ2')

sc.pp.neighbors(scscWTg2.sc_matrix, n_neighbors=40, n_pcs=20, use_rep='X_pca_harmony')  #
sc.tl.umap(scscWTg2.sc_matrix, min_dist=0.4)
cmap = clr.LinearSegmentedColormap.from_list('3', ['#E6E6E6', '#54B7DF', '#006497', '#022234'], N=256)
sc.pl.umap(scscWTg2.sc_matrix, color=['cellType', 'dpt_pseudotime', 'irootCluster'], show=False, cmap=cmap,  save='GermAfterDiff')

newOrder = ['1', '6', '3', '4', '0', '5', '2']
scscWTg2.cell_type_identification(name='leiden', order=newOrder, markerGeneFile='markersGerm2') #)
scscWTg2.cell_type_identification(clusterMethod='leiden', name='GermZF', order=newOrder, markerGeneFile='markersCompareZFgerm') #)

# scscWTg.marker_genes(cluster='leiden', name='WTGerm')

# ####################Male####################
scscM = copy.deepcopy(scsc)
scscM.sc_matrix = scsc.pick_specific_cells('Group', ['WTM', 'KOM'])
cells = scscWT.sc_matrix[scscWT.sc_matrix.obs['cellType'] == 'Somatic'].obs_names.to_list()
cells.extend(scscM.sc_matrix[scscM.sc_matrix.obs['Group'] == 'KOM'].obs_names.to_list())
scscM.sc_matrix = scscM.pick_specific_cells_by_name(cells)
scscM.filtering_normalization_scaling()

scscM.pca_and_umap(neighbors=20.0, pcs=20.0, res=0.4, dist=0.5, name='M', size=20)
newOrder = ['10', '7', '3', '9', '6', '0', '4', '8', '2', '1', '5']
scscM.cell_type_identification(name='M', order=newOrder, size=20,  markerGeneFile='markersMale')

# scscM.marker_genes(cluster='leiden', name='M')
"""
# differential expression
# scscM.DEgenes = []
scscM.sc_matrix.obs['leidenGroup'] = (scscM.sc_matrix.obs['leiden'].astype(str) + scscM.sc_matrix.obs['Group'].astype(str)).astype('category')  #().astype('category')
for i in ['0', '1', '2', '5']: # no sub-sampling
# for i in ['0', '1', '2', '4']:      # sub-sampling
    scscM.differential_expression_genes(i+'WTM', i+'KOM', cluster='leidenGroup', name='M_%s_WTvsKO' % i, thresholdpval=0.1, thresholdFC=1.5)
scscM.sc_matrix = scscM.pick_specific_cells('leiden', ['0', '1', '2', '5'])
scscM.FCheatmap(genesViolin, name='male_FC_heatmap')


####################Female#################
scscF = copy.deepcopy(scsc)
scscF.sc_matrix = scsc.pick_specific_cells('Group', ['WTF', 'KOF'])
cells = scscWT.sc_matrix[scscWT.sc_matrix.obs['cellType'] == 'Somatic'].obs_names.to_list()
cells.extend(scscF.sc_matrix[scscF.sc_matrix.obs['Group'] == 'KOF'].obs_names.to_list())
scscF.sc_matrix = scscF.pick_specific_cells_by_name(cells)

scscF.filtering_normalization_scaling()
# scscF.sc_matrix = scscF.subsampling_by_group('KOF', 'WTF')
# print(pd.crosstab(index=scscF.sc_matrix.obs['Group'], columns='count'))

scscF.pca_and_umap(neighbors=20.0, pcs=20.0, res=0.8, dist=0.5, name='F', size=20)
newOrder = ['15', '13', '11', '16', '7', '14', '8', '9',  '0', '1', '4', '5', '10', '2', '3', '17', '12', '6']
scscF.cell_type_identification(name='F', order=newOrder, size=20, markerGeneFile='markersFemale')
print(pd.crosstab(index=scscF.sc_matrix.obs['Group'], columns='count'))
# scscF.marker_genes(cluster='leiden', name='F')


# # differential expression
scscF.sc_matrix.obs['leidenGroup'] = (scscF.sc_matrix.obs['leiden'].astype(str) + scscF.sc_matrix.obs['Group'].astype(str)).astype('category')
print(pd.crosstab(index=scscF.sc_matrix.obs['leidenGroup'], columns='count'))

# for i in ['0', '1', '2']:             # sub-sampling
    scscF.differential_expression_genes(i+'WTF', i+'KOF', cluster='leidenGroup', name='F_%s_WTvsKO' % i, thresholdpval=0.1, thresholdFC=1.5)

# scscF.sc_matrix = scscF.pick_specific_cells('leiden', ['0', '1', '2', '5'])
# # scscF.sc_matrix = scscF.pick_specific_cells('leidenGroup', ['5WTF', '5KOF', '2WTF', '2KOF', '0WTF', '0KOF', '1WTF', '1KOF'])

