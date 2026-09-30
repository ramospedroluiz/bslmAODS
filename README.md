# bslmAODS

R code for **Bayesian Supervised Learning Methods: Regularization and Predictive Uncertainty**, by Pedro Luiz Ramos. Accepted for publication in *Annals of Data Science*.

The illustration predicts concrete compressive strength using the same 1,030 observations and five composition-grouped folds as the article. All repository files are in this directory. No manuscript, installed packages or large fitted objects are included.

## Run

Open a terminal in this directory. R 4.6.0 is the reference version.

```sh
Rscript --vanilla setup.R
Rscript --vanilla run.R
```

The first command restores the recorded package versions from `renv.lock`. It needs internet access and installs packages in an R user-data library **outside this repository**. Compilation tools may be needed if a matching binary is unavailable. Set `BSLM_LIBRARY` to an existing library path to use a different location.

The second command regenerates the score tables, bootstrap intervals and calibration figure from the saved row-level forecasts. It checks the numerical results against the reference files. It **does not refit the models**. For this quick route, `ggplot2` is the only non-base R package required, so an existing installation can also be used.

```sh
Rscript --vanilla run.R check   # Also verify the data and grouped folds
Rscript --vanilla run.R test    # Also run predictive and short sampler tests
Rscript --vanilla run.R fit     # Refit all models, score, summarize and test
Rscript --vanilla run.R score   # Rescore completed local fits without refitting
```

The full fit is substantially slower and produces several GB of local fitted objects and checkpoints. It requests one computational thread. Do not run two full fits in the same directory at once. BART resumes compatible checkpoints; changes to the data, fitting scripts or recorded package environment prevent their reuse. If the maximum sampling budget fails the diagnostic thresholds, scoring stops instead of reporting that fit as validated.

## Files

| File | Purpose |
| --- | --- |
| `run.R` | Main entry point; defaults to reproducing saved results |
| `setup.R`, `renv.lock` | Package restoration |
| `config.R` | Settings, training-only standardization and flat file paths |
| `prepare.R` | Reconstruct and verify the grouped folds from the original spreadsheet |
| `models.R`, `fit-models.R` | Eight non-tree primary fits and predictive distributions |
| `fit-bart.R` | Final 1,000-tree BART specification, diagnostics and continuation |
| `sensitivity.R` | BMA prior sensitivity; GAM covariance sensitivity is scored in `evaluate.R` |
| `evaluate.R`, `summarize.R` | Row-level forecasts, scores, bootstrap summaries and figure |
| `tests.R`, `test-bart.R`, `verify.R` | Statistical implementation and reproduction checks |
| `Concrete_Data.xls`, `concrete.csv`, `DATA.md` | Original data, prepared data with fixed fold IDs, attribution and column definitions |
| `reference-*` | Frozen results and diagnostic/environment records from the article |
| `PROTOCOL.md` | Analysis decisions and chronological BART amendments |

Generated tables and plots start with `output-`, for example `output-summary.csv` and `output-calibration.pdf`. Fits start with `fit-` and BART checkpoints with `chain-`. These files stay in the same directory and are ignored by Git. Reference files are never overwritten. The R library and temporary installation files are stored outside the repository.

## Methods and interpretation

The primary comparison contains OLS (classical reference), Bayesian ridge, Bayesian lasso, horseshoe, Bayesian model averaging, GAM (REML), BART, GP (empirical Bayes) and a neural network (evidence/Gauss-Newton approximation). GAM, GP and the neural network do not integrate all hyperparameter uncertainty; they should not be described as exact full Bayesian fits.

BART uses **1,000 trees**, four chains per fold, 8,000 warm-up iterations, and 16,000 retained draws per chain initially. The same chains may be continued to 32,000 and 64,000 retained draws. The diagnostic gate is maximum R-hat <= 1.01 and minimum bulk and tail ESS >= 400 over the monitored summaries. The machine-readable name `BART-1000` corresponds to `BART` in the article.

The saved forecasts also retain the BMA/GAM sensitivity analyses and earlier BART fits explicitly labelled `diagnostic_only`. The latter are excluded from the primary results. The simplified fitting route implements the final specification, not the historical unsuccessful BART pilots. Those amendments, including their timing relative to inspection of the initial results, are documented in `PROTOCOL.md`.

Bootstrap intervals resample 428 compositions using **fixed out-of-fold forecasts**; models are not refitted in each resample. The comparison is illustrative, not a claim of universal model superiority. Historical zero elapsed-time readings for OLS/BMA are below the timer's resolution, not literally zero computational cost. Saved-object sizes are not peak memory measurements.

## Reproducibility

The `reference-*` files are unchanged exports from the accepted article's computational material. In particular, the source hashes in `reference-bart-execution.json` identify the **original** scripts, not the renamed and reorganized scripts here. New fits record their own source/data/package identifiers. Original and final-BART session-information records were identical; one copy is included.

The flat workflow was checked by regenerating the primary score table, composition-bootstrap intervals and paired comparisons, verifying data/fold isolation, checking predictive mixture identities, and testing native BART output and serialized chain continuation. Full fitting was not repeated merely to reorganize this repository. The short BART test is an implementation test, not a convergence assessment of a new analysis. Package restoration was previously tested for the article; it is not a guarantee that every future platform or package server will behave identically. Full refits may differ numerically across platforms or library builds.

## Citation and licenses

Please cite the accompanying article when referring to this analysis:

Pedro Luiz Ramos. *Bayesian Supervised Learning Methods: Regularization and Predictive Uncertainty*. Annals of Data Science, accepted for publication.

Repository: <https://github.com/ramospedroluiz/bslmAODS>. `CITATION.cff` supplies the software citation. No article DOI or pagination is asserted here before it is available.

No license has been assigned to the repository's code at the author's request. This does not change the UCI dataset's **CC BY 4.0** license or the licenses of third-party R packages. See `DATA.md` for the dataset citation and attribution.
