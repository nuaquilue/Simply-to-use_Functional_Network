## Simply-to-use_Functional_Network
Data and R code to cluster species into functional groups, and compute stand-level functional diversity and vulnerability 

### Input files

##### species.txt
List of species
**code**: unique four letters code for species
**species**: species name

##### functional.trait.txt
List of nine functional trait values for 75 tree and shrub species
**code**: unique four letters code for species
**seed.mass**: Seed mass (g·1000-1 seeds)
**wood.dens**: Wood density (g·cm-3)
**leaf.mass.area**: Leaf mass area (g·cm-2)
**shade.tol**: Index of shade tolerance (1 - intolerant to 5 - tolerant)
**drought.tol**: Index of drought tolerance (1 - intolerant to 5 - tolerant)  
**waterlog.tol**: Index of water logging tolerance (1 - intolerant to 5 - tolerant)
**phylogen**: Phylogenetic division (A - Angiosperm, G - Gymnosperm)
**dispersion**: Dispersal vector (A - Animal, H - Human, U - Unassisted, W - Wind)
**ass.mycorrhiza**: Mycorrhiza association (0 - no, 1 - yes)

##### spp.rel.abund.stand.txt
Species relative abundance per forest stand or tree community
**code**: unique id for stands
**ABBA** - **TSCA**: relative abundance of species code as in **species**

##### scores.txt
Raw scores (-3 - negativelly influenced to 3 - positivelly influenced) for a list of species and natural disturbances

##### disturbances.txt
List of disturbances and relative uncertainity and future relevance
**disturbance**: unique disturbance name
**future.relev**: future relevance (1 - low, 2 - moderate, 3 - high, 4 - very high)
**uncertainty**: uncertainty (0 to 1)
