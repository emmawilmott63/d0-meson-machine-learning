# Directory Layout

```
macros/      C++ and ROOT macros
scripts/     run_generation.sh
data/        (created)
plots/       (created)
logs/        (created)
tmva/        (created)
  plots/     (created)
```

# Commands

```
make       # build programs (same as make all)
make all   # build programs 
make clean # delete programs
```
### Run Individual Programs

Simulate p+p collision events in PYTHIA.
```
make generate EVENTS=1000 JOBS=3   
# or 
./scripts/run_generation.sh 1000 EVENTS 3
```
Output: `data/d0_inc_5_100.root`

Create a ROOT TTree with K-Pi pair candidates.
```
make candidates 
# or 
./buildtree_jetbins data/d0_inc_5_100.root
```
Output: `data/candidates_inc_5_100.root`

Create a ROOT TTree with only opposite-charge K-Pi pair candidates.
```
make candidates-os
# or 
./buildtree_jetbins_os data/d0_inc_5_100.root
```
Output: `data/candidates_os_inc_5_100.root`

Plot the TTree candidate distribution for each recorded feature.
```
make features-plots
./features_jetbins data/candidates_inc_5_100.root
```

```
make features-plots-os
./features_jetbins_os data/candidates_os_inc_5_100.root
```

Confirm TMVA (Toolkit for Multivariate Analysis) compatibility before running.
```
make check-tmva
```

Run TMVA machine learning with BDT (boosted decision trees).
```
make train-bdt-jetbins
# or 
root -l -b -q 'train_bdt_jetbins.C("data/candidates_inc_5_100.root","tmva")'
# or (opposite-charge candidates only)
root -l -b -q 'train_bdt_jetbins.C("data/candidates_os_inc_5_100.root","tmva")'
```