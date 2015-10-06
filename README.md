# TM_density
Python script to calculate transmembrane density of particles, residues and charges around transmembrane proteins.
TM profiles are calculated with respect to the local bilayer normal.
To print this help menu: ```python TM_density.py --help```

```
DESCRIPTION
 
 This script plots the density profile along the z axis of groups of particles present
 around transmembrane proteins. Density profiles are broken down by the size of the
 TM clusters.

 A file containing the charged particles can also be supplied to calculate the density
 of charges.
 
 density calculation
 -------------------
 Density is calculated based on a cylinder centered on the protein cluster center of 
 geometry. All the particles present in this cylinder are binned into cylinder slices
 based on their distance (along z) to the center of the lipid bilayer. The cylinder
 radius is controlled via --slices_radius and the slices thickness via --slices_thick.

 (a) particles selection
   You can specify the particles for which to plot the density by supplying a file
   via the --particles option. Each line of this file should follow the following
   format (without quotation marks):
    -> 'group,label,colour,MDAnalysis selection string'
 
   where 'group' is used to normalise the densities of several particles. By default
   each particle type is normalised with respect to itself and the following densities
   are used:
    -> protein,protein,#262626,protein
    -> CHOL,CHOL,#bdbdbd,resname CHOL and name ROH
    -> POPC,POPC,#41ab5d,resname POPC and name PO4
    -> POPE,POPE,#6a51a3,resname POPE and name PO4
    -> POPS,POPS,#cc4c02,resname POPS and name PO4
    -> water,water,#1d91c0,resname W
    -> Na+,Na+,#7bccc4,name NA+
    -> Cl-,Cl-,#fa9fb5,name CL-
  
   Note that the only acceptable protein selection is the whole "protein" string and
   that it should be labeled as "protein" and in a group of its own also labeled
   "protein". If you're interested in particular sub-selections of proteins, see (b)
   below.
   
 (b) protein details: residue types
   Density profiles can be broken down further for peptides to show the density
   distribution amongst different residue types. You can group the residues as you
   wish by supplying a file via the --residues options. Each line of this file
   should follow the following format:
   -> 'label,colour,resname1,resname2,...'
  
   By default the following residue types are used:
    -> basic,#253494,ARG,LYS
    -> polar,#a1d99b,SER,THR,ASN,GLN,HIS		
    -> negative,#99d8c9,ASP, GLU
    -> hydrophobic,#993404,VAL,ILE,LEU,MET,PRO,CYS,PHE,TYR,TRP
    -> backbone,#969696,ALA,GLY

   WARNING: you need to have defined a particle type called "protein" in (a) in order
   to show residue details.
   If you do not want to show residue details, just use: '--residues no'.

 (c) charge density
   You can specify which particles to take into account for the calculation of the total
   charge density by supplying a file via the --charges option. Each line of this file
   should follow the format (without quotation marks):
    -> 'group_name,colour,charge_name,charge,MDAnalysis selection string for charge'

   Note that the MDAnalysis selection string should not contain any commas.

   The total charge for each group will be plotted on the charge density profile. The
   group colour must be specified for each charge.
  
   By default charges are defined using the Martini 2.1 settings:
    -> solvent,#52A3CC,Na+,1,name NA+
    -> solvent,#52A3CC,CL-,-1,name CL-
    -> lipids,#b2182b,-1,phosphate,name PO4
    -> lipids,#b2182b,1,amine_choline,name NH3 or name NC3
    -> protein,#053061,Lys,1,resname LYS and name SC2
    -> protein,#053061,Arg,1,resname ARG and name SC2
    -> protein,#053061,Asp,-1,resname ASP and name SC1
    -> protein,#053061,Glu,-1,resname GLU and name SC1
  
   Another default set of charges can be used by specifiying --charges 2.2P :
    -> solvent,#52A3CC,Na+,1,name NA+
    -> solvent,#52A3CC,CL-,-1,name CL-
    -> solvent,#52A3CC,WP,0.46,name WP
    -> solvent,#52A3CC,WM,-0.46,name WM
    -> lipids,#b2182b,-1,phosphate,name PO4
    -> lipids,#b2182b,1,amine_choline,name NH3 or name NC3    
    -> protein,#053061,Lys,1,resname LYS and name SCP
    -> protein,#053061,Arg,1,resname ARG and name SCP
    -> protein,#053061,Asp,-1,resname ASP and name SCN
    -> protein,#053061,Glu,-1,resname GLU and name SCN
    -> protein,#053061,Asn_p,0.46,resname ASN and name SCP
    -> protein,#053061,Asn_n,-0.46,resname ASN and name SCN
    -> protein,#053061,Gln_p,0.46,resname GLN and name SCP
    -> protein,#053061,Gln_n,-0.46,resname GLN and name SCN
    -> protein,#053061,Thr_p,0.31,resname THR and name SCP
    -> protein,#053061,Thr_n,-0.31,resname THR and name SCN
    -> protein,#053061,Ser_p,0.4,resname SER and name SCP
    -> protein,#053061,Ser_n,-0.4,resname SER and name SCN

   By default proteins termini are considered uncapped and +1 and -1 are added to the
   first and last backbone beads ("name BB") respectively. If this is not what you want
   just use the option --capped.
   
   Note that for systems containing different protein species automatic addition of
   charges at the termini has not yet been implemented: you should use --capped and
   provide a file via --charges specifying all charged particles in the system, including
   termini beads.

   If you do not want to calculate charge density, just use: '--charges no'

 (d) colour definition
   Colours can be specified using single letter code (rgbcmykw), hex code  or the name of
   a colour map (see the matplotlib website for a list of the available colour maps).
   In case a colour map is used, its name must be specified as the only colour.

 detection of transmembrane protein clusters
 -------------------------------------------
 Two clustering algorithms can be used to identify protein clusters.

 (a) Connectivity based (relies on networkX module):
   A protein is considered in a cluster if it is within a distance less than --nx_cutoff
   from another protein. This means that a single protein can act as a connector between
   two otherwise disconnected protein clusters.
   This algorithm can be ran using either the minimum distante between proteins (default, 
   --algorithm 'min') or the distance between their center of geometry (--algorithm 'cog').
   The 'min' option scales as the square of the number of proteins and can thus be very
   slow for large systems.

 (b) Density based (relies on the sklearn module and its implementation of DBSCAN):
   A protein is considered in a cluster if is surrounded by at least --db_neighbours other
   proteins within a radius of --db_radius.
   This density based approach is usually less suited to the detection of protein
   clusters but as a general rule the more compact the clusters, the smaller --db_radius
   the higher --db_neighbours can be - for details on this algorithm see its online
   documentation.
   This algorithm is selected by setting the --algorithm option to 'density'.

 The identified protein clusters are considered to be transmembrane only if the closest
 lipid headgroup neighbours to the cluster particles are all within the same leaflet.
 In addition to the sizes identified, size groups can be defined - see note 7.


REQUIREMENTS

 The following python modules are needed :
  - MDAnalysis
  - matplotlib
  - numpy
  - scipy
  - networkX (if option --algorithm is set to 'min' or 'cog')
  - sklearn (if option --algorithm is set to 'density')


NOTES

 1. The density is calculated with respect to the z axis, not the bilayer normal. So the
    more your system deforms the noiser the less meaningful the results get.

 2. Identification of the bilayer leaflets can be controlled via 3 options:
    (a) beads
     By default, the particles taken into account to define leaflet are:e
     -> name PO4 or name PO3 or name B1A
   
     Note that only lipids which contain one of the beads mentioned in the selection string
     will be taken into account. If you wish to specify your own selection string (e.g. to
     choose different beads or add a bead not in the default list in order to take into
     account a particular lipid specie) you can do so by supplying a file via the --beads
     option. This file should contain a single line that can be passed as the argument
     to MDAnalysis selectAtoms() routine and should not contain any quotation marks, e.g.:
      -> name PO4 or name PO3 or name B1A or name AM1
        
    (b) leaflet finding method
     By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
     the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
     routine.
     This optimisation process can take time in large systems and you can specify your own
     cutoff value to skip this step. For instance to use a 15 Angstrom cutoff value:
      -> '--leaflet 15'
    
     In very large systems (more then ~50,000 phospholipids) LeafletFinder (or rather the
     networkX module that it relies on) can fail. To  avoid this you can choose not to use
     this routine by specifying:
      -> '--leaflet large'
     In this case lipids whose headgroups z value is above the average lipids z value will
     be considered to make up the upper leaflet and those whose headgroups z value is below
     the average will be considered to be in the lower leaflet.
     This means that the bilayer should be as flat as possible in the 1st frame of the xtc
     file supplied in order to get a meaningful outcome. 

	 NOTE: By default the gro file is only used as a topology file and the 1st frame of the
	 xtc is used to identify leaflets. If you wish to use the gro file instead, for instance
	 in the case that the 1st frame of the xtc is not flat, you need to specify the --use_gro
	 flag: be warned that this might take a few minutes longer on large systems.

    (c) flipflopping lipids
     In case lipids flipflop during the trajectory, a file listing them can be supplied
     with the --flipflops option. Each line of this file should follow the format:
      -> 'resname,resid,starting_leaflet,z_bead'
     where starting_leaflet is either 'upper' or 'lower' - e.g. 'POPC,145,lower,PO4'. The
     z_bead is used to track the position of the lipid.
     If flipflopping lipids are not specified they may add significant noise to results as
     they prevent the correct identification of TM clusters.
     Note that, even when specified, flipflopping lipids will be taken into account when
     calculating densities and charges.   

 3. Proteins are detected automatically but you can specify an input file to define your
    own selection with the --proteins option.
    In this case the supplied file should contain on each line a protein selection string
    that can be passed as the argument of the MDAnalysis selectAtoms() routine - for 
    instance 'bynum 1:344'.
   
 4. The densities are calculated for each TM cluster size identified but can also be
    binned into size groups.
    The size groups are defined by supplying a file with --groups, whose lines all
    follow the format:
     -> 'lower_size,upper_size'

    Size groups definition should follow the following rules:
     -to specify an open ended group use 'max', e.g. '3,max'
     -groups should be ordered by increasing size and their boundaries should not overlap
     -boundaries are inclusive so you can specify one size groups with 'size,size,colour'
     -any cluster size not covered will be labeled as 'other'
 
 5. There are 3 possible options to determine the local normal to the bilayer. These are
    controlled with the flags --normal and --normal_d:
    (a) 'z': the bilayer is assumed flat in the x,y plane and the z axis is taken to be the
     normal. Best for systems without significant curvature and local deformations. In this
     case the --normal_d flag is ignored.

    (b) 'cog': in this case neighbourhing particles to current cluster of interest are
     identified in the lower and upper leaflet. The local normal is then considered to be the
     vector going from the cog of the lower ones to the cog of the upper ones. In each leaflet,
     neighbouring particles are the particles selected by --beads which are within --normal_d
     Angstrom of the cog of the protein cluster of interest.

    (c) 'svd': in this case neighbourhing particles to current cluster of interest are
     identified in the lower and upper leaflet as in (b) above. The normal of the best fitting
     plane to these particles is obtained by performing a singular value decomposition of their
     coordinate matrix.

 
USAGE

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc]
-o			: name of output folder
-b			: beginning time (ns) (the bilayer must exist by then!)
-e			: ending time (ns)	
-t 		10	: process every t-frames
--particles		: definition of particles, see 'DESCRIPTION'
--types			: definition of residue groups, see 'DESCRIPTION'
--charges	2.1	: definition of charged particles, see 'DESCRIPTION' 
--capped		: assumes protein termini are capped
 
Density profile options
-----------------------------------------------------
--range		40 	: distance spanned on either side of the bilayer center
--slices_thick	0.5 	: z thickness of the slices (Angstrom)
--slices_radius	30 	: radius of the slices (Angstrom)
--normal	svd	: local normal to bilayer ('z', 'cog' or 'svd'), see note 5
--normal_d	50	: distance of points to take into account for local normal, see note 5

Lipids identification  
-----------------------------------------------------
--beads			: leaflet identification technique, see note 2(a)
--flipflops		: input file with flipflopping lipids, see note 2(c)
--leaflets	optimise: leaflet identification technique, see note 2(b)
--use_gro		: use gro file instead of xtc, see note 2(b)

Protein clusters identification
-----------------------------------------------------
--groups		: cluster groups definition file, see note 4
--proteins		: protein selection file, (optional, see note 3)
--algorithm	min	: 'cog','min' or 'density', see 'DESCRIPTION'
--nx_cutoff 	6	: networkX cutoff distance for protein-protein contact (Angstrom)
--db_radius 	20	: DBSCAN search radius (Angstrom)
--db_neighbours	3	: DBSCAN minimum number of neighbours within a circle of radius --db_radius	

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
```
