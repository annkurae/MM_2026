import argparse
import gc
import logging
import sys
from pathlib import Path
from typing import Optional, Callable

import anndata as ad
from omegaconf import OmegaConf

project_root = Path(__file__).resolve().parent.parent 
sys.path.append(str(project_root / 'src'))

from step4_infer_annotation import (
    infer_missing_major_population,
    infer_tme_population,
    infer_pc_population
)
from io_utils import generate_path_in_output_dir
from logging_utils import set_file_logger
from train_scvi_model import load_pp_adata_after_norm_and_hvg


def infer_annotation_of_all_cells(config, save_path: Path, on_adata: Optional[ad.AnnData] = None):
    if on_adata is None:
        adata = load_pp_adata_after_norm_and_hvg(config)
    else:
        adata = on_adata
    drop_diseases = ('In_vitro', 'Ex_vivo')
    logging.info(f"dropping {drop_diseases} from {config.annotation.Disease} column")
    adata = adata[adata.obs[config.annotation.Disease].apply(lambda x: x not in drop_diseases)].copy()

    logging.info("inferring major population, PC vs TME, using generic SCVI model")
    infer_missing_major_population(config, adata)

    logging.info("inferring TME sub population, using generic SCVI model")
    infer_tme_population(config, adata)

    logging.info("inferring PC sub population, using generic SCVI model")
    infer_pc_population(config, adata)

    if save_path is not None:
        logging.info(f"saving AnnData with inferred annotation to file - {save_path}")
        adata.write(save_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='infer annotation of all cells',
        description='create h5ad file from pp data, with full annotation')

    default_config = str(project_root / 'configs' / 'config.yaml')
    parser.add_argument('--config', help='a path to an valid config file', default=default_config)

    args = parser.parse_args()

    conf = OmegaConf.load(args.config)

    logging_file_path = Path(conf.outputs.output_dir, conf.outputs.logging_file_name)
    set_file_logger(logging_file_path, prefix="infer")

    output_path = generate_path_in_output_dir(conf, 'adata_with_scvi_full_annot_pred.h5ad',
                                              add_version=True, add_date_timestamp=True)
    infer_annotation_of_all_cells(config=conf, save_path=output_path)
