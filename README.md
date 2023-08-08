# TNBC
The project combines spatial statistics and ABM to look at CSF effects on the remission-relapse dynamics of TNBC. The spatial stats are used to see if 
there is any bias in the spatial distribution of proliferating cells with respect to their distance to the stroma. 
There is in vitro and indirect in vivo evidence that stroma (CAF) produces paracrine factors that enhance proliferation, hence in stroma's proximity, a higher rate of proliferation would 
outcompete the death to due chemotherapy and results in relapse without the acquisition of hard-wired resistance. 
We use a proliferation marker that stains cells in the S phase of the  mitosis only (BrdU) to identify cells that are undergoing active cell division and use
spatial statistics to asses their spatial relationship with respect to the stroma. 
Specifically, the mouse TNBC xenografts are collected at different time points along the chemotherapy treatment, and the IHC samples are converted into 
point patterns using Aiforia, an AI platform, and R (PxToCoordTNBC.R). We use these point patterns to find spatial summary statistics metrics, mainly using 
the spatstat package in R (RanalysisWholeTissue.R and RanalysisWholeTissueRand.R)
We then extract the point pattern for the stroma and populate the ABM grid with, while each point in the grid inherits the information about how many grid points away
the grid point is from the nearest stroma grid (allDist/allDist4G). This information will determine if the tumor cell residing there will benefit from stroma sheltering or not.
The grid is then populated with cells at all locations, and it's let run the simulation (ExampleGrid.java) until it equilibrates to about 20% of the grid point being empty (based on the empty
space noticed in the actual IHC tissues due to micronecrosis and migration).
