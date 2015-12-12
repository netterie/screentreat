# screentreat

## Version guide for Annals paper
- Click on "version_guide_Annals.csv" to view the list of model versions stored in the GitHub repository and some distinguishing details.
- The first column indicates the corresponding Results table in the Annals paper
- The results contained in the Results tables in the Annals paper are taken from the output file "examples/[insert-model-version-name]/output/cuminc_mrr_newtable.csv"

## How to run and understand the models
- Clone this repository to your local machine
- Select a model to run locally (breast_ER-HER2_6 is the main model of the Annals paper)
- Follow instructions in [run_README.txt](https://github.com/netterie/screentreat/blob/master/run_README.txt) to run the model locally. 
    - In Step 1, you will additionally need to edit the _rootdir_ object in the "directory setup" section in the wrapper file to point to the directory in which you have cloned the screentreat repository.
    - In the user_options file of the model, you may initially want to decrease the _nsim_ object from 100 to a smaller number, e.g. 2, until you're confident everything is working.
- Once you have successfully run the model, step through the code to understand the inputs and processes.

## How to choose which model to run
1. See version_guide_AnnalsPaper.xlsx to select a model corresponding to results in the Annals paper
