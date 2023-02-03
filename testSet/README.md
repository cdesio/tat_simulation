# RBE Test Set
The aim of this test set is to provide examples which can be run to ensure setup is working. The scripts run 3 separate simulations which can be compared to a reference set. The reference set was run for 1000 primaries.

Prior to running this test set the programs should be build following the main README.

Make a results folder in testSet. All commands should be run from ./testSet/results.

## simulation and clustering
This will run three simulations and run the clustering algorithm on all three. This takes approx 30 mins locally, but only runs 100 primaries for each simulation. On blue pebble 1000 primaries are run.

Firstly plotChromatinFibreSection.py must be run to create the geometry files.

Check the following in plotChromatinFibreSection.py before running:
- volumeWidth = 300 
- numRows = 4 
- numColumns = 4

run (python3):
```
python ../../geometryFiles/plotChromatinFibreSection.py
```

Ensure that the results file only contains the 2 geometry files before running the simulations.

If running locally run the following to run the simulations (takes approx. 30 mins):
```
chmod 755 ../runAllSims.sh (this will make the script executable)

../runAllSims.sh
```

If running on blue pebble run the following to run the simulations:
```
sbatch ../runSim1.sh

sbatch ../runSim2.sh

sbatch ../runSim3.sh
```

To check the progress of simulations on blue pebble run the following to see if they are still running:
```
squeue -u $USER$
```
This should take approximately 1 hour on blue pebble. 

## Compare to expected results
When the simulations are finished three sets of csv files will be created (simX.csv and DSBclusterSize_simX.csv). run plotResults.py to compare the simulation results to previous results.
```
python ../plotResults.py
```



