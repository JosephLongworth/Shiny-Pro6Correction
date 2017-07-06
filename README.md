# Shiny-Pro6Correction
Shiny script creating a web application for the correction of SILAC ratios generated from Maxquant for proline 6 incorporation.


The easiest way to run the web application is through R as shown below. Alternatively the Server.R and Ui.R can be downloaded and run locally or uploaded to a server as described at http://shiny.rstudio.com/.

```R
# 1. If not installed first install the following packages from CRAN.
install.packages("shiny")
install.packages("stringr")

# 2. Load Shiny.
library(shiny)

# 3. Run the application from GitHub.
runGitHub("Shiny-Pro6Correction", "JosephLongworth")


