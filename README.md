# Covariates_Analysis

This repository is linked to the publication: 

[La Fisca et. al., "Biases in BCI Experiments: Do we really need to balance stimulus properties across categories?", Frontiers in Computational Neuroscience, Volume 16, Pages 153-172, 22 November 2022](https://www.frontiersin.org/articles/10.3389/fncom.2022.900571/full)

All descriptions and computations of covariates related to image visual features can be found in the [utils](https://github.com/numediart/Covariates_Analysis/blob/main/utils/Readme.md) folder.

## Requirements
- [FieldTrip](https://github.com/fieldtrip/fieldtrip.git)
- [LIMO-for-FieldTrip](https://github.com/LucaLaFisca/LIMO-for-FieldTrip.git)
- [limo_tools](https://github.com/LIMO-EEG-Toolbox/limo_tools.git)

If you want to run the code as it is, you should download the corresponding dataset from [Zenodo](https://zenodo.org/records/7298746#.Y2kKIXbMK3A)

## Tutorial
### 1. Configuration of the study
Adapt the config.json file to your data.

a) Information about the study:
```
	"model_names"		 	: "model_cat,model_psycho,model_image,model_psycho_image",
	"model_names_naive" 	: "model_psycho,model_image,model_psycho_image",
	"save_choice" 			: true,
	"SOURCE_ANALYSIS" 		: false,
	"task_name" 			: "semanticPriming",
	"datatype"				: "bdf",
	"comment_datatype"		: "bdf, edf or mat",
```

b) Information about the paths:
```
	"PATH_TO_ROOT" 			: "D:/__EEG-data/BIDS_files",
	"PATH_TO_FIELDTRIP" 	: "D:/FieldTrip",
	"PATH_TO_LIMO" 			: "D:/_LIMO/limo_tools",	
	"PATH_TO_FT2LIMO" 		: "D:/_LIMO/LIMO-for-FieldTrip",
	"PATH_TO_SOURCE" 		: "D:/__EEG-data/BIDS_source",
```

c) Information about the confounders:
```
	"PATH_TO_COV_ANALYSIS" 	: "C:/Users/luca-/OneDrive - UMONS/_PhD/__GitHub/Covariates_Analysis",
	"PATH_TO_IMAGES" 		: "D:/_ARC-images/ARC_tache_Final/300x300",
	"PATH_TO_SIMILARITY" 	: "utils/visual_similarity.csv",
	"PATH_TO_ITEMS" 		: "utils/items_description.csv",
  "cov_description" 		: "utils/continuous_variables.csv",
	"comment_cov_descript." : "give the path to the csv file containing the values of the studied confounders",	
```

d) Information about the trials:
```
	"trial_function" 		: "SenSem_trialfun_trial",
	"comment_trial_function": "replace with the desired trial function (cf. FieldTrip doc)",
	"trialdef_eventtype"	: "STATUS",
	"trialdef_eventvalue" 	: [10, 11, 12],
	"trialdef_prestim"		: 0.5,
	"trialdef_poststim"		: 1,
	"trial_start" 			: -200,
	"trial_end" 			: 500,
	"comment_trial" 		: "in ms",
```

e) Information for the source analysis:
```
	"atlas"					: "template/atlas.mat",
	"headmodel"				: "template/vol.mat",
	"elec"					: "template/elec.mat",
	"channels"				: "[1:64]"
```

### 2. Selection of the confounders to study
The function ```var_selection(config)```, that can be found in [utils](https://github.com/numediart/Covariates_Analysis/blob/main/utils/Readme.md) folder, performs a correlation analysis of the different confounders defined in the csv file provided in "cov_description".
