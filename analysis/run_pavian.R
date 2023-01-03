# Pavian isn't a super common package. Here is how to download it if you don't have it
if (!require(remotes)) { install.packages("remotes") }
if (!require(pavian)) { remotes::install_github("fbreitwieser/pavian") }

pavian::runApp(port=5000)
