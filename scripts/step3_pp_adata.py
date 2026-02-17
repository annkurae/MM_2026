import argparse
import logging
import os
import sys
from pathlib import Path
import gc

import anndata as ad
from omegaconf import OmegaConf, DictConfig

# sys.path.append(os.getcwd())
# this assumes your script is in a folder like 'scripts/' or 'notebooks/'
# at the top level of project
project_root = Path(__file__).resolve().parent.parent 
sys.path.append(str(project_root / 'src'))
# sys.path.append(str(project_root / 'scripts'))

from data_loading.utils import load_dataframe_from_file, merge_labels_to_adata
#from 2_load_sc_data_to_anndata import load_sc_data_to_anndata
from data_loading.load_data_to_anndata import load_sc_data_to_anndata as load_sc_data_to_anndata
from logging_utils import set_file_logger
from pre_processing.black_list import drop_blacklist_genes_by_pattern
from pre_processing.scanpy_pp import pp_drop_genes, pp_drop_cells

from io_utils import generate_path_in_output_dir


def load_annotations(adata: ad.AnnData, config: DictConfig):
    annotation_df = load_dataframe_from_file(Path(config.annotation.annotations_file_name))
    merged_adata = merge_labels_to_adata(adata, annotation_df, col_in_adata_to_merge_by="index",
                                         cols_in_labels_df_to_merge_by=config.annotation.cell_id_columns_name,
                                         cols_to_validate_not_empty=[], merge_suffixes=(None, "_annotation"))
    adata.obs = merged_adata.obs
    if "Healthy" in set(adata.obs.Disease) and "Healthy " in set(adata.obs.Disease):
        adata.obs.Disease = adata.obs.Disease.apply(lambda x: "Healthy" if x == "Healthy " else x)


def pre_process(config: DictConfig):
    adata_path_from_config = generate_path_in_output_dir(config, config.outputs.loaded_adata_file_name,
                                                         add_version=True)
    if adata_path_from_config.exists():
        logging.info(f"reading raw AnnData from {adata_path_from_config}")
        adata = ad.read_h5ad(adata_path_from_config)
    else:
        adata = load_sc_data_to_anndata(config)

    if config.annotation.annotations_file_name != "None":
        logging.info(f"loading annotations from {config.annotation.annotations_file_name}")
        load_annotations(adata, config)

    if config.annotation.exclude_blood:
            initial_cells = adata.shape[0]
            tissue_col = 'Tissue'  
            if tissue_col not in adata.obs.columns:
                raise ValueError(
                    f"'{tissue_col}' column not found in adata.obs. "
                    f"Available columns: {list(adata.obs.columns)}"
                )
            
            tissue_dist_before = adata.obs[tissue_col].value_counts()
            logging.info(f"Tissue distribution before exclusion:\n{tissue_dist_before}")
            
            excluded_tissues = ['Blood', 'In_vitro', 'Ex_vivo']
            adata = adata[~adata.obs[tissue_col].isin(excluded_tissues)].copy()
            cells_removed = initial_cells - adata.shape[0]
            logging.info(
                f"Excluded cells from {', '.join(excluded_tissues)}: {cells_removed} cells removed "
                f"({adata.shape[0]} cells remaining)"
            )
            gc.collect()

    logging.info(f"dropping genes from blacklist")
    adata = drop_blacklist_genes_by_pattern(adata, config.pp.blacklist_genes_pattern_list)
    pp_drop_genes(adata, min_num_cells_per_gene=config.pp.min_num_cells_per_gene)
    pp_drop_cells(adata, min_num_genes_per_cell=config.pp.min_num_genes_per_cell,
                  min_num_counts_per_cell=config.pp.min_num_counts_per_cell,
                  mitochondrial_patterns=config.pp.mitochondrial_gene_patterns,
                  max_pct_mt_genes_pre_cell=config.pp.max_pct_mt_genes_pre_cell)
    logging.info(f"dropping mitochondrial genes")
    adata = drop_blacklist_genes_by_pattern(adata, config.pp.mitochondrial_gene_patterns)

    processed_adata_path = generate_path_in_output_dir(config, config.outputs.processed_adata_file_name,
                                                       add_version=True)
    if processed_adata_path is not None:
        adata.write(processed_adata_path)
        logging.info(f"saving processed AnnData to file - {processed_adata_path}")
    return adata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='AnnData Loading',
        description='loads scRNA data to AnnData and save it to h5ad file')

    default_config = str(project_root / 'configs' / 'config.yaml')
    parser.add_argument('--config', help='a path to an valid config file', default=default_config)

    args = parser.parse_args()

    conf = OmegaConf.load(args.config)

    logging_file_path = generate_path_in_output_dir(conf, conf.outputs.logging_file_name)
    set_file_logger(logging_file_path, prefix="pp")

    pre_process(config=conf)
