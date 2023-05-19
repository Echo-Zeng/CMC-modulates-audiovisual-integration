# CMC-modulates-audiovisual-integration


## Experimental Data
*data_total_exp1.mat* contains behavioral data from Experiment 1, *data_total_exp2.mat* contains behavioral data from Experiment 2. Both files share the following column tags:
- tms_type: -1 for no TMS stimulation, 1 for TMS stimulation P4, 2 for TMS stimulation Cz
- block: 0 for auditory only condition, 1 for physical visual stimulation, 2 for visual mental imagery
- v_type: -1 for no visual stimulus, 1 for risky visual stimulus, 2 for round visual stimulus
- loc_v: NaN for no visual stimulus condition, other values refers to the horizontal locations of visual stimuli in degrees
- a_type: 0 for the sound "kiki", 1 for the sound "bouba"
- loc_v: values refers to the horizontal locations of auditory stimuli in degrees
- discrepancy: NaN for no visual stimulus condition, other values refers to spatial disparities between visual and auditory stimuli
- cmc: crossmodal correspondence, 1 for congruent crossmodal correspondence, 0 for incongruent crossmodal correspondence
- response: reported perceived auditory location
- rt: reaction time

## Bayesian Causal Inference (BCI) Model fitting
*bciBasicContinuous.m* and *bciMultiCondition.m* are main functions to perform BCI model fitting. These two functions are customized based on the [BCIT by Ladan Shams](https://github.com/multisensoryperceptionlab/BCIT) and the [codes shared by Ferrari and Noppeney (2021)](https://doi.org/10.1371/journal.pbio.3001465).
