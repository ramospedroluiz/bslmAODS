# Data

Yeh, I-Cheng (1998). **Concrete Compressive Strength** [Data set]. UCI Machine Learning Repository. <https://doi.org/10.24432/C5PK67>.

Source: <https://archive.ics.uci.edu/dataset/165/concrete+compressive+strength>.
License: [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).
Downloaded for the analysis on 9 September 2026. The dataset is redistributed with attribution; UCI and the data creator do not endorse this analysis.

`Concrete_Data.xls` is the original spreadsheet. `concrete.csv` preserves its 1,030 observations, renames the nine original columns and appends composition, row and fold identifiers. All 25 exact repeated rows are retained. No response values are changed.

| Column | Meaning |
| --- | --- |
| `cement`, `slag`, `flyash` | Cement, blast-furnace slag and fly ash (kg/m3) |
| `water`, `superplasticizer` | Water and superplasticizer (kg/m3) |
| `coarse`, `fine` | Coarse and fine aggregate (kg/m3) |
| `age` | Age (days); log-transformed during preprocessing |
| `strength` | Concrete compressive strength (MPa), the response |
| `recipe_id` | Identifier for an identical seven-ingredient composition |
| `row_id` | Original row number |
| `fold` | Fixed grouped cross-validation fold (1 to 5) |

All observations with the same composition, including different ages, stay together. The 428 compositions are assigned to five folds with 206 rows each without using the response. Standardization is estimated separately within each training fold. Composition grouping does not identify production batches or near-duplicate mixtures.

`prepare.R` reconstructs the prepared table from the spreadsheet and verifies it without overwriting the distributed data. SHA-256 of the original spreadsheet:

```text
710076c66b9ca3f8050e7942f3dcbdbe04013534daeb0077ffd3079a52d8e0c4
```

In `reference-predictions.csv`, `y` is the observed strength; `mean` and `variance` are predictive moments in MPa and MPa squared; `qNNN` is the predictive quantile at probability NNN/1000. `nlpd`, `crps` and `pit` denote negative log predictive density, continuous ranked probability score and probability integral transform. `coverNN` and `widthNN` describe the central NN-percent predictive interval. `role` distinguishes primary, sensitivity and diagnostic-only forecasts. `recipe_id`, `row_id` and `fold` link each prediction to `concrete.csv`.
