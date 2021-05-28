MATLAB files as referenced in "Analyzing the differences in Olfactory Bulb spiking with ortho- and retronasal stimulation"
Folder structure:
-Analytic = MATLAB scripts for synaptic input statistic equations and linear filter model. Files:
     -check_SynMnVar_MC2MCt.m = Script to compare Monte Carlo ORN input statistics to analytical ODE apprx
     -drvLinFiltShift.m = Driver script of LinFiltShift.m to fit inputs (time-varying syn) to MonteCarlo of ouput stats
-OBsc = Olfactory Bulb single compartment biophysical model. Files:
     -[Orth/Retr]Bt.m = Driver function to iterate multiple trials of the OB Single Compartment model (obsc[O/R]t.m) for Ortho/Retro stimulus.
     NOTE: This script is written to only generate a single trial because the model is computationally expensive.
     The statistics in d[Orth/Retr]_ct.mat files are averaged over 50,000 trials.
     -d[Orth/Retr]_ct.mat = Workspace of OB network model first and second-order statistics averaged over 50,000 trials
     -expDat_orEB_Rat1.mat = Workspace of in vivo rat data first and second-order statistics (single rat, multiple cells and cell pairs)
     -show_OrRetr.m = Script to compare OB network model MC stats with Orth/Retr in vivo data (uses .mat files)
