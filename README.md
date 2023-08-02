# GEOS-Chem Kinetic PreProcessor (KPP) Files
### Species and equations for the full chemistry mechanisms in KPP v2.3.3 (GEOS-Chem Classic v13.3.3)

**custom_eqn.txt**: includes the Automated Model Reduction (AMORE) v1.1 isoprene oxidation mechanism 

**fullchem_eqn.txt**: includes the default isoprene oxidation mechanism (BASE)

# Implementing the AMORE mechanism in KPP
1. Navigate to your KPP directory (in Linux)
```
cd ~/GCClassic/src/GEOS-Chem/KPP
```
2. Download the "custom.eqn" file from GitHub
```
wget https://github.com/benjaminyang93/amore_v1.1/blob/main/custom.eqn
```
4. Overwrite the existing "custom" mechanism file
```
mv custom.eqn custom/custom.eqn
```
6. Run KPP to create new chemical solver Fortran-90 files
```
./build_mechanism.sh custom
```
For more https://kpp.readthedocs.io/en/stable/

5. Recompile model (cmake ../CodeDir -DRUNDIR=.. -DCUSTOMMECH=y)

### For more information about how to use KPP and GEOS-Chem, please see the following guides:
https://kpp.readthedocs.io/en/stable/index.html
https://geos-chem.readthedocs.io/en/stable/
