READ ME	

Updated and final code versions after review will be uploaded to a public repository.

Code to take and process the raw CITES data product can be found at Morton et al. (2022) Curr Bio.

Obj1_Model_and_outputs 
- Fits the data Traded_Traits_Data.csv and processes the outputs to reproduce all analyses ascociated with Obj1 in the main paper.
- NOTE phylogenetic models take some time to run on a standard PC.

Obj2_Model_Fitting
- This fits CITES_fitting_data.csv (note this is the cleaned data so can be fitted straight to the model without the cleaning steps in the code as it is already cleaned).
- NOTE phylogenetic models take some time to run on a standard PC.

Obj2_Phylosignal
- Used both datasets and calculates the phylosignal and reproduces Figure 2.

Obj2_Plotting_and_Interpretation
- Uses the model output of Obj2_Model_Fitting and CITES_fitting_data.csv
- reproduces figure 3 and 4 and writes out all relevant summaries for model interpretation.