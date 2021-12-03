# ATLAS
The original ATLAS program (http://www.cbr.washington.edu/analysis/apps/atlas) provides the ability to
analyze mark-recapture data that used active-tag technology, and to adjust the observed survival estimates
for potential tag-failure that may be misidentified as a mortality of the marked animal. The cbrATLAS R
package was created to allow this analysis be conducted in R. In addition, cbrATLAS incorporates an
update of the original version to include the R package failCompare, increasing the number of model fitting
options of tag-life data to predict tag-life of the active tags used in the study. The failCompare package also
provides goodness-of-fits and rankings of the various models, to assist selection of the most appropriate
model of tag-life, along with plotting the expected tag-life curve against observed tag failures.

cbrATLAS is currently capable of analyzing single release data histories (with and without censoring
events) that may or may not require correction for active-tag failure. Future iterations are planned to
incorporate paired release , Virtual-Paired Release (ViPRe), Virtual Release Dead Fish Correction
(ViRDCt) study designs.

Please see the vignette for examples of a single release mark-recapture analysis, with and without correction for potential tag failure.
