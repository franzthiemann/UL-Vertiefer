# UL-Vertiefer
Containes Scripts developed during the Internship in the Tonner-Zech Group at the University of Leipzig

## amsresults.py
Can be executed in a folder that that was used to perform pEDA or NOCV calculation using the pEDA script in a slurm job.
The script will parse the pEDA command in all slurm filee, evaluate the pEDA calculations using pEDA_eval. In addition to that, NOCVs are rendered using VMD and the results are written into a summary file, that containes all pEDA and NOCV results with their respective path.
The script requires no arguments, however, the paths for the VMD template file, vmdrc file and summary file need to be set prior to usage.


## interpolateMixed.py
Used to generate NEB-path by mixing the geosidic interpolation for the molecule with an IDDP interpolation for the surface.

### Usage
```
interpolateMixed.py [start] [end] [Images] [SurfaceCutoff] <seperate>
```
| start: | POSCAR file of initial molecule |
| --- | --- |
| end:  |          POSCAR file of final molecule |
| Images:  |       Number of NEB Images to generate |
| SurfaceCutoff: | Atom Number of the first Atom that belongs to the molecule |
| seperate:  |     Pass 'True' when convergence is difficult to ommit the molecule from the surface interpolation |

requires the geosidic_interpolate package to be installed (https://github.com/virtualzx-nad/geodesic-interpolate)
