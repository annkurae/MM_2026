from pathlib import Path

import pandas as pd

if __name__ == '__main__':
    data_version = '2025-06-03'

    plates_metadata_excel_path = Path(f"/home/projects/amit/annaku/data/mm_2023/Blueprint_MM_Plates_{data_version}.xlsx")
    sheet_names = ["MM3_MARS", "MM3_SPID", "Blood_MARS", "Blood_SPID"]

    data_version = plates_metadata_excel_path.stem.split("_")[-1]

    for sheet_name in sheet_names:
        tissue = "Blood" if "Blood" in sheet_name else "BM"
        method = "MARS" if "MARS" in sheet_name else "SPID"
        csv_path = Path(plates_metadata_excel_path.parent,
                        f"Blueprint_MM_{tissue}_{method}_Plates_{data_version}.csv")
        df = pd.read_excel(plates_metadata_excel_path, sheet_name=sheet_name)
        df["Method"] = method

        df.to_csv(csv_path)
