#!/bin/bash
cd /Users/plavskin/Documents/yeast\ stuff/Mutational\ Effect\ Modeling/plavskin_code/pipeline_v2
matlab -nodisplay -nosplash -nodesktop -nojvm -r 'matlab_wait('"30"");exit" | tee ${output}/logs/PIFluor_out.txt