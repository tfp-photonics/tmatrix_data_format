## ----eval=TRUE, message=FALSE-------------------------------------------------
library(glue)
library(dplyr)
library(purrr)

# read in the template
template <- glue_collapse(readLines('../smarties/batch/_template.m'), sep = "\n")

parameters <- expand.grid(a = seq(10, 100, by=5), c = seq(10, 100, by=5),
                          material = c("Au", "Ag", "Si"), 
                          medium = c("vacuum", "water")) |> 
  filter(a != c) |> # use Mie for this
  mutate(shape = ifelse(a > c, "oblate", "prolate"),
         n = ifelse(medium == "water", 1.33, 1.0),
         source = case_when(material == "Au" ~ "Raschke et al 10.1103/PhysRevB.86.235147", 
                            material == "Ag" ~ "Raschke et al 10.1103/PhysRevB.91.235137",
                            material == "Si" ~ "Aspnes et al 10.1103/PhysRevB.27.985",
                            .default = "Unknown material!"))

parameters$step <- 1
nrow(parameters)


## ----eval=TRUE----------------------------------------------------------------
parameters <- expand.grid(a = seq(20, 50, by=10), c = seq(20, 50, by=10),
                          material = c("Au", "Ag"), 
                          medium = c("water")) |> 
  filter(a != c) |> # use Mie for this
  mutate(shape = ifelse(a > c, "oblate", "prolate"),
         n = ifelse(medium == "water", 1.33, 1.0),
         source = case_when(material == "Au" ~ "Raschke et al 10.1103/PhysRevB.86.235147", 
                            material == "Ag" ~ "Raschke et al 10.1103/PhysRevB.91.235137",
                            material == "Si" ~ "Aspnes et al 10.1103/PhysRevB.27.985",
                            .default = "Unknown material!"))
parameters$step <- 5


## ----echo=FALSE---------------------------------------------------------------
gt::gt(parameters[,c("a","c","shape","material","medium")])


## -----------------------------------------------------------------------------
write_script <- function(a, c, material, medium, shape, n, step, source){
 script <- glue(template)   
 cat(script, file = glue('../smarties/batch/run_{material}_{medium}_{a}_{c}.m'))
 cat(glue("run_{material}_{medium}_{a}_{c}\n\n"), 
     file = "../smarties/batch/batch.m", append = TRUE)
}
cat("%% Running all the files below\n",file = "../smarties/batch/batch.m", append = FALSE)
pwalk(rowwise(parameters), write_script)

