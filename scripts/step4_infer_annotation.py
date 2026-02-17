import argparse
import logging
from datetime import date
from pathlib import Path
import sys
from typing import Optional, Callable

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from omegaconf import OmegaConf
from sklearn.neighbors import KNeighborsClassifier

project_root = Path(__file__).resolve().parent.parent 
sys.path.append(str(project_root / 'src'))

from io_utils import generate_path_in_output_dir
from logging_utils import set_file_logger
from sc_classification.utils import add_number_of_patients_in_neighborhood, count_number_of_annotation_in_neighborhood
from train_scvi_model import load_pp_adata_after_norm_and_hvg, generate_model_name


def _get_latest_model_like(config, model_name_with_no_timestamp: str) -> Path:
    #logging.info(f"Looking for models matching pattern: {model_name_with_no_timestamp}")
    all_fitting_model_paths = list(
        Path(config.outputs.output_dir).glob(f"{model_name_with_no_timestamp}_20*"))  # old format
    all_fitting_model_paths += list(Path(config.outputs.output_dir).glob(f"{model_name_with_no_timestamp}_ts_20*"))
    if len(all_fitting_model_paths) == 0:
        raise ValueError(
            f"no model in output dir: {config.outputs.output_dir}, fitting pattern: {model_name_with_no_timestamp}_20*")
    
    #logging.info(f"Found {len(all_fitting_model_paths)} matching models")

    all_fitting_suffixes = [str(model_path).split("_")[-1] for model_path in all_fitting_model_paths]
    all_fitting_dates = []
    for suffix in all_fitting_suffixes:
        try:
            all_fitting_dates.append(date.fromisoformat(suffix))
        except ValueError:
            all_fitting_dates.append(date.min)
    model_path = all_fitting_model_paths[np.argmax(all_fitting_dates)]
    #logging.info(f"Selected most recent model: {model_path}")
    return model_path


def _load_model_and_compute_latent_and_neighborhood(adata, config, model_path, scvi_latent_key, neighborhood_key,
                                                    compute_umap=True, compute_leiden=False,
                                                    add_normalized_expression=False,
                                                    transform_batch=None, transform_key="scvi_expr",
                                                    leiden_resolution: Optional[float] = None):
    logging.info(f"loading model: {model_path}")
    model = scvi.model.SCVI.load(str(model_path), adata=adata)
    adata.obsm[scvi_latent_key] = model.get_latent_representation()
    logging.info(f"computing neighborhood to key: {neighborhood_key}")
    sc.pp.neighbors(adata, use_rep=scvi_latent_key,
                    n_neighbors=config.umap_settings.knn_k, key_added=neighborhood_key)
    if compute_leiden:
        res = config.annotation_prediction.global_manifold.leiden_resolution if leiden_resolution is None else leiden_resolution
        logging.info(f"computing leiden clustering with resolution: {res}")
        sc.tl.leiden(adata, resolution=res, neighbors_key=neighborhood_key)
    if compute_umap:
        logging.info(f"computing umap with neighborhood key: {neighborhood_key}")
        sc.tl.umap(adata, min_dist=config.umap_settings.umap_min_dist, neighbors_key=neighborhood_key)
    if add_normalized_expression:
        logging.info(f"computing normalized expression with batch: {transform_batch}, to layer key: {transform_key}")
        adata.layers[transform_key] = model.get_normalized_expression(adata, n_samples=10,
                                                                      return_mean=True, transform_batch=transform_batch)


def _update_predicted_annotations_in_adata(adata, config, annotation_col, y_pred_series):
    adata.obs[annotation_col].update(y_pred_series)
    predicted_col_name = config.annotation_prediction.indication_col_name
    adata.obs[predicted_col_name] = adata.obs.index.isin(y_pred_series.index)


def _infer_adata_annotation_on_scvi(adata: ad.AnnData, config, annotation_column, knn_neighborhood_key, scvi_model_path,
                                    single_cell_emb_model, scvi_latent_key,
                                    adata_cell_condition: Optional[Callable] = None):
    logging.info(f"scvi settup for anndata of shape: {adata.shape}")
    scvi.model.SCVI.setup_anndata(adata, layer=config.scvi_settings.counts_layer_name,
                                  batch_key=config.scvi_settings.batch_key)
    if scvi_latent_key not in adata.obsm.keys():
        logging.info(f"computing latent representation to {scvi_latent_key}")
        _load_model_and_compute_latent_and_neighborhood(adata, config, scvi_model_path, scvi_latent_key,
                                                        knn_neighborhood_key, compute_leiden=True)
    else:
        logging.info("using precomputed embedding")

    inference_adata = adata.copy()
    if adata_cell_condition is not None:
        mask = inference_adata.obs.apply(adata_cell_condition, axis=1)
        logging.info(f"conditioning: dropping cells before inference, before : {len(mask)}, after: {sum(mask)}")
        inference_adata = inference_adata[mask]

    logging.info(f"infering values of column: {annotation_column}")
    adata_with_annot = inference_adata[~inference_adata.obs[annotation_column].isna()]
    adata_cells_to_predict = inference_adata[inference_adata.obs[annotation_column].isna()]

    X_train = adata_with_annot.obsm[scvi_latent_key]
    y_train = adata_with_annot.obs[annotation_column]
    X_pred = adata_cells_to_predict.obsm[scvi_latent_key]
    logging.info(f"triaing annotation model to predict: {np.unique(y_train)}")
    single_cell_emb_model.fit(X_train, y_train)
    y_pred = single_cell_emb_model.predict(X_pred)

    y_pred_series = pd.Series(y_pred, index=adata_cells_to_predict.obs_names)
    _update_predicted_annotations_in_adata(adata, config, annotation_column, y_pred_series)


def infer_missing_major_population(config, adata: ad.AnnData):
    model_path = _get_latest_model_like(config, generate_model_name(config))
    logging.info(f"Using model for major population inference: {model_path}")
    pc_model = KNeighborsClassifier(n_neighbors=config.annotation_prediction.pc_vs_tme.k)

    _infer_adata_annotation_on_scvi(adata, config,
                                    scvi_model_path=model_path,
                                    single_cell_emb_model=pc_model,
                                    annotation_column=config.annotation.major_cell_type_column,
                                    knn_neighborhood_key=None if config.scvi_settings.neighborhood_key == "None" else config.scvi_settings.neighborhood_key,
                                    scvi_latent_key=config.scvi_settings.scvi_latent_key)

    return adata


def _infer_annotation_on_sub_population(config, adata, adata_cell_condition: Callable):
    scvi_model_path = _get_latest_model_like(config, generate_model_name(config))
    logging.info(f"Using model for major population inference: {scvi_model_path}")
    pc_model = KNeighborsClassifier(n_neighbors=config.annotation_prediction.within_tme.k)
    _infer_adata_annotation_on_scvi(adata, config,
                                    scvi_model_path=scvi_model_path,
                                    single_cell_emb_model=pc_model,
                                    adata_cell_condition=adata_cell_condition,
                                    annotation_column=config.annotation.cell_type_column,
                                    knn_neighborhood_key=None if config.scvi_settings.neighborhood_key == "None" else config.scvi_settings.neighborhood_key,
                                    scvi_latent_key=config.scvi_settings.scvi_latent_key)


def infer_tme_population(config, adata):
    def cd45_condition_in_row(row):
        return row[config.annotation.major_cell_type_column] == "CD45"

    _infer_annotation_on_sub_population(config, adata, cd45_condition_in_row)


def infer_pc_population(config, adata):
    def pc_condition_in_row(row):
        return row[config.annotation.major_cell_type_column] == "PC"

    _infer_annotation_on_sub_population(config, adata, pc_condition_in_row)


def _filter_pc_with_sure_annotation(adata: ad.AnnData, config) -> ad.AnnData:
    count_number_of_annotation_in_neighborhood(adata, config.annotation.major_cell_type_column, "PC")
    count_number_of_annotation_in_neighborhood(adata, config.annotation.major_cell_type_column, "CD45")
    pcs_with_majority_tme_neighbors = (adata.obs[config.annotation.major_cell_type_column] == "PC") & (
            adata.obs['count_of_PC_in_neighborhood'] < adata.obs['count_of_CD45_in_neighborhood'])

    cluster_annot = {}
    for c, group in adata.obs[["leiden", config.annotation.major_cell_type_column]].groupby("leiden"):
        cluster_annot[c] = group.value_counts().index[0][1]
    pcs_in_tme_clusters = (adata.obs[config.annotation.major_cell_type_column] == "PC") & \
                          (adata.obs['leiden'].apply(lambda x: cluster_annot[x] == "CD45"))

    pc_with_tme_environment = pcs_with_majority_tme_neighbors | pcs_in_tme_clusters
    adata.obs["pc_with_tme_environment"] = pc_with_tme_environment
    only_pc_adata = adata[(adata.obs[config.annotation.major_cell_type_column] == "PC") & (~ pc_with_tme_environment)]

    return only_pc_adata


def process_pc_population(config, adata, save_only_pc_path: Optional[Path] = None, return_only_pc=False):
    description = f"_only_{'_'.join(['Malignant', 'Healthy'])}"
    pc_scvi_model_path = _get_latest_model_like(config, generate_model_name(config, extra_description=description))
    logging.info(f"model path: {pc_scvi_model_path}")

    _load_model_and_compute_latent_and_neighborhood(adata, config,
                                                    model_path=pc_scvi_model_path,
                                                    scvi_latent_key=config.scvi_settings.only_pc_scvi_latent_key,
                                                    neighborhood_key=config.scvi_settings.only_pc_neighborhood_key,
                                                    transform_batch=config.scvi_settings.counts_imputation_batch,
                                                    transform_key=config.scvi_settings.pc_expression_layer,
                                                    compute_umap=False,
                                                    add_normalized_expression=False)

    if (save_only_pc_path is not None) or return_only_pc:
        logging.info(f"computing neighbors umap and leiden for only pc")
        only_pc_adata = _filter_pc_with_sure_annotation(adata, config)
        sc.pp.neighbors(only_pc_adata, use_rep=config.scvi_settings.only_pc_scvi_latent_key,
                        n_neighbors=config.umap_settings.knn_k,
                        key_added=config.scvi_settings.only_pc_neighborhood_key)
        sc.tl.umap(only_pc_adata, min_dist=config.umap_settings.umap_min_dist,
                   neighbors_key=config.scvi_settings.only_pc_neighborhood_key)
        sc.tl.leiden(only_pc_adata, resolution=config.annotation_prediction.within_pc.leiden_resolution)
        add_number_of_patients_in_neighborhood(only_pc_adata, patient_col=config.annotation.patient_id_column_name,
                                               distances_obsp_key=f"{config.scvi_settings.only_pc_neighborhood_key}_distances")
        if save_only_pc_path is not None:
            only_pc_adata.write(save_only_pc_path)
            logging.info(f"saving only pc AnnData with clustering to file - {save_only_pc_path}")
        if return_only_pc:
            return only_pc_adata


def infer_annotation_of_unlabeled_cells(config, on_adata: Optional[ad.AnnData] = None, save_path=None):
    if on_adata is None:
        adata = load_pp_adata_after_norm_and_hvg(config)
    else:
        adata = on_adata
    drop_diseases = ('In_vitro', 'Ex_vivo')
    logging.info(f"dropping {drop_diseases} from {config.annotation.Disease} column")
    adata = adata[adata.obs[config.annotation.Disease].apply(lambda x: x not in drop_diseases)].copy()


    logging.info("inferring major population, PC vs TME, using generic SCVI model")
    infer_missing_major_population(config, adata)

    # return None
    logging.info("inferring TME sub population, using generic SCVI model")
    infer_tme_population(config, adata)

    logging.info("process PC sub population, will not infer but do all processing before, using PC SCVI model")
    if save_path is not None:
        save_path = Path(save_path)
        pc_save_path = save_path.with_stem(
            save_path.stem + "_only_pc")  # TODO move the '_only_pc' between the name and timestamp
    else:
        pc_save_path = None

    process_pc_population(config, adata, save_only_pc_path=pc_save_path)

    if save_path is not None:
        logging.info(f"saving AnnData with inferred annotation to file - {save_path}")
        adata.write(save_path)

    return adata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='infer annotation of cells',
        description='create h5ad file from pp data, with partial annotation')

    default_config = str(project_root / 'configs' / 'config.yaml')
    parser.add_argument('--config', help='a path to an valid config file', default=default_config)

    args = parser.parse_args()

    conf = OmegaConf.load(args.config)

    logging_file_path = Path(conf.outputs.output_dir, conf.outputs.logging_file_name)
    set_file_logger(logging_file_path, prefix="infer")

    adata_full_annotations_path = generate_path_in_output_dir(conf, conf.outputs.inferred_missing_annotation_file_name,
                                                              add_version=True, add_date_timestamp=True)
    infer_annotation_of_unlabeled_cells(config=conf, save_path=adata_full_annotations_path)
