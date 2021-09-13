**Codes to reproduce the results in the paper "Modeling animal movement with directional persistence and attractive points"**

Models estimates can be obtained running the Julia script *"PosteriorEstimates.jl"*. Inside the script, a value must be given to the objects *"IndModel"*, *"NAME"* and *"DIR"*.


 * *"DIR"* is the directory where all the relevant codes and folders are saved;
 * *"IndModel"* is a vector with 1 element, that has value 1 if we want to estimate the STAP-HMM, value 2 for the BRW-HMM, and 3 for the CRW-HMM. 
 * The results, depending on the value of *"IndModel"*, are saved in an R object called *"NAME_STAP.Rdata"*, *"NAME_BRW.Rdata"* or *"NAME_CRW.Rdata"* in the directory *"DIR/out/"*.

To automatically load the exact versions of the Julia packages, the file *"DIR/Juliapackagesversions/Manifest.toml"* should be loaded and activated before running the code. Inside the file, *"XXX"*  on line 32 must be replace with *"DIR/BayesianAnimalMovementModels.jl"*

All figures and tables can be reproduced running the script *"Plots.R"*. Inside the script, a value must be given to the objects *"PLOT_DIRDATA"*, *"NAME"* and *"PLOT_DIRPLOT"*, *"MOD_STAP_NAME"*, *"MOD_BRW_NAME"* and *"MOD_CRW_NAME"*.

 * *"PLOT_DIRDATA"* is the directory where the results of the script *"PosteriorEstimates.jl"* are saved (*"DIR/out/"*),
 * *"PLOT_DIRPLOT"* is the directory where we want to save the figures,
 * *"MOD_STAP_NAME"*, *"MOD_BRW_NAME"* and *"MOD_CRW_NAME"* are the names of the R objects that contain the posterior samples of the respective model.

Tables are printed on the R console. To have a clearer output, the file *"PlotsAndTables.R"* should be run in R using the command *"source("DIR/PlotsAndTables.R")"*

Notice that there are inconsistencies between the names of the parameters in the paper and in the scripts: *\mu* and *\nu* in the paper are, respectively, *\mu_0* and *\psi* in the scripts.

The folder *"BayesianAnimalMovementModels.jl"* contains the Julia package needed to estimate the models.

The folder *"Data"* contains the dataset saved in a *".Rdata"* format (R workspace); coordinates of several dogs and sheep are present.
The one used in the paper are *"Long_Bindi_Canis"* (Longitude) and  "Lat_Bindi_Canis" (Latitude).

The ".csv" data can be downloaded at "https://www.datarepository.movebank.org/handle/10255/move.395".
