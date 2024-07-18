# utilities for generating the website

library(knitr)

# for some reason purl() comments out the code when run here
# extract the code from each qmd file as standalone script
# purl("website/export_julia.qmd",  "website/export_julia.jl")
# purl("website/export_python.qmd", "website/export_python.py")
# purl("website/export_R.qmd",      "website/export_R.R")
# purl("website/export_matlab.qmd", "website/export_matlab.m")
# 
# purl("website/smarties_advanced.qmd", "smarties/example_advanced.m")
# purl("website/smarties_simple.qmd", "smarties/example_simple.m")

# copy scripts to utilities for main branch
file.copy("website/export_julia.jl", "utilities/export_julia.jl")
file.copy("website/export_python.py", "utilities/export_python.py")
file.copy("website/export_R.R", "utilities/export_R.R")
file.copy("website/export_matlab.m", "utilities/export_matlab.m")

file.copy("website/smarties_simple.m", "smarties/smarties_simple.m")
file.copy("website/smarties_advanced.m", "smarties/smarties_advanced.m")
file.copy("website/smarties_batch.m", "smarties/smarties_batch.m")

