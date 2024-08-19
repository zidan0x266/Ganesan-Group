# utanalysis
## Ions Transport analysis tools for polyelectrolytes

The code uses MDAnalysis as the basic analyzer, currently it works with Gromacs format trajectory, it should also work with Lammps format trajectory without any problem but without guaranty. **Python3 needed**

The code have four functions now:
- analyze the association relationship of a given ion to the particles of interests
- analyze the association behavior of ions near interface
- analyze the dynamic properties of the ion transport system, including intermittent autocorrelation functions C(t) and continuous autocorrelation function S(t)
- analyze the hopping event types of ion transport

Coming functions:
- local RDF in individual domain when morphology is taken into account
- display the multiprocessing tracking bar by using **tqdm** package
- introduce the chunksize to improve the performance of multiprocessing, since the pool.map() function assigns the job sequencially while it is always ture that the first few jobs are more time-consuming

## Usage
After obtaining the trajectory file from either atomistic __AT__ or coarse-grained __CG__ simulation (for atomistic simulation, we recommend that the user firstly transfer the trajectory to __CG__ with the mapping scheme as the center of mass of the functional groups), at least two groups of atoms need to be defined prior to anything else, namely the anion group **aP** and cation group **aN**, then the analysis should work automatically as using the following commands. 
For those system with interface, an additional group need to be defined to track the interface, namely **aE**.

- Association analysis
```
python3 utanalysis.py --top topol.tpr --trj traj.xtc --begin 100 --end 200 --func homopoly --nt 20 --cutoff 6.3 --counterN 8 --anion PF --cation IM
```
- Hopping events analysis
```
python3 utanalysis.py --func hopping --h5 pils.h5 --h5tag htt --counterN 8
```
- Association time-autocorrelation function analysis
```
python3 utanalysis.py --func aacf --h5 pils.h5 --h5tag htt --trestart 10 --typeacf ct --nt 20
```

__OR__
```
python3 utanalysis.py @params
```
User cound find params sample file in examples
### Install with pipenv

You can very simply install all required dependencies by using `pipenv`.

```
pipenv shell
```

will install all required dependencies and move you to the shell where you can run the command
from the **Usage** section

### Indications
- top - topology file
- trj - trajectory file
- begin - frame to begin with, default = 0
- end - frame to end with
- func - copoly: random/block copolymer; homopoly: homo polymer, default = copoly
- nt - number of processors
- trestart - time shifting window for calculating aacf
- h5 - input file for hopping and aacf
- h5tag - dataset name in .h5 file
- typeacf - type of the acf to be calculated
- anion - anion type
- cation - cation type
- beadiny - bead in the second domain to determine the interface
