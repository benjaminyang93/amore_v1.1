## GEOS-Chem Kinetic PreProcessor (KPP) Chemical Mechanism Files
### Species and equations in KPP v2.3.3 (GEOS-Chem Classic v13.3.3)

* **custom_eqn.txt**: includes the Automated Model Reduction (AMORE) v1.1 isoprene oxidation mechanism 

* **fullchem_eqn.txt**: includes the default isoprene oxidation mechanism (BASE)

## Implementing the AMORE mechanism in KPP and GEOS-Chem 
1. Navigate to your KPP directory (in Linux)
```
cd <path>/GCClassic/src/GEOS-Chem/KPP
```
2. Download the "custom.eqn" file from GitHub
```
wget https://raw.githubusercontent.com/benjaminyang93/amore_v1.1/main/custom.eqn
```
3. Overwrite the existing "custom" mechanism file
```
mv custom.eqn custom/custom.eqn
```
4. Run KPP to create new chemical solver Fortran-90 files
```
./build_mechanism.sh custom
```
5. Navigate to your GEOS-Chem build directory and recompile
```
cmake ../CodeDir -DRUNDIR=.. -DCUSTOMMECH=y; make -j install
```
6. Run the model faster than with the BASE mechanism!  

### For more information on how to use KPP and GEOS-Chem, please see the following:

https://kpp.readthedocs.io/en/stable/index.html

https://geos-chem.readthedocs.io/en/stable
