# Gromacs Job Production/Analysis

This code aims at providing tools to easily generate gromacs simulation protocol and analyse results with parts based on already existing codes.

## Getting Started

The following instructions will get you a copy of the project up and running on your local machine.

### Prerequisites

The whole code is written in Python3 so you will need python3 package on your machine.

In addition to Python3 standard library, the following packages are needed for Production.py:

```
numpy
packmol
vmd
```
And for Analysis.py:

```
numpy
pandas
matplotlib
scipy
mdanalysis
fatslim
```



#### On Debian

```
sudo (or su -c) apt-get install python3-pip
```

And then for any system having [pip3](https://pip.pypa.io/en/latest/):

```
pip3 install --upgrade numpy pandas matplotlib scipy mdanalysis fatslim
```

with vmd and packmol available with the following links:

[vmd](https://www-s.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)

[packmol](http://m3g.iqm.unicamp.br/packmol/userguide.shtml)

### Installing

To get the code on your machine simply use the following command:

```
git clone https://github.com/fben94/GJB
```

## How to use

### Production

To create simulation protocols for gromacs a .csv file must be written similarly to what is shown in the default Parameters.csv .
A first block describes the name of the project, the paths to gromacs (on local and remote machine), packmol and vmd and finally the path to Gromacs default files.
Gromacs default files contain topology files for MARTINI, default molecular dynamic parameters (.mdp), base membranes/solvent/substrate and topology plus MD parameters for defo, substrate and wall.

After this first block jobs are represented in two parts.
A first part describes the system to simulate, job related parameters and gromacs mdrun options (see Parameters.csv with optional sections tagged with #).
The second part starting after the section **PROTOCOL** describes the simulation protocol. The first rwo after the protocol section describes the way the system is generated. It is either initialised with the parameters set in the first part using the **INIT** section or obtained from another simulation using the **INPUT** section.
After that, each two rows correspond to a simulation step with its associated parameters and their set values. The first column sets the name of the step and the next columns define the parameters and their values. The format is repeated for every step.

The overall format system + simulation steps is repeated for every job and separated using empty rows. To stop reading the file before its end a **END_OF_PROJECT** section can be set before other jobs.

This next part describes the different sections for a job:

#### JOBID

This section describes the name, number and packmol seed for the current job.

#### NODE

With this section you can set the name of the node number of node and processor per node to used.
On this same row, an optional parameter `SEQUENTIAL` can be set to have simulations steps run as different jobs (useful when runs are time limited).

#### MDRUN_OPT

This section pass parameters to mdrun for all the simulation step.

#### SAMPLE (for INIT)

The system molecules and their associated number are described here.

#### TYPE (for INIT)

The type of the simulation is set here and can be *TRILAYER*, *BILAYER*, *MONOLAYER* and *SOLVENT*. Each type has different optional parameters for the geometry section in addition to required parameters.

#### GEOMETRY (for INIT)

This section sets the geometry of the box. Required parameters are the box dimensions `LX`, `LY`, `LZ` in Angstrom. Optional parameters depend on system type with the previous section. 

For this section there is no parameter order.

#### DEFO (to include a pore in the membrane)

This section aims at setting the artficial molecule (DEFO) parameters.

- The `Height` which can be set to *follow*, *bilayer* or *box*. The former creates defo molecules for each membrane (bilayer and monolayer for example), the second only for the bilayer and the latter along the whole box in the z-direction.
- The distance between layers of artificial particles (DEF) is set using `DzDefo` in Angstrom.
- The number of DEF in a layer is set with `DpL` (DEF per layer) with 0 meaning only the DEF at the center.
- For packmol sample generation, the radius beyond which the lipid positions can be initialised is set using the `Radius` parameter in nanometre.
- The topology version of the DEF particles is set using the `Version` parameter.
- To set constraint between DEF you can set the parameter `Constraints` to *bonds*, *angles* or *bond&angles*.
- The force between in-plane, out-of-plane bonds and angles can be set using `FbondN`, `FbondP`, `FangleN`, and `FangleP` resp. with the value in kJ/mol/nm.
- The type of bond and angle potentials for gromacs is set with `FbondType` and `FangleType` resp.
- The last parameter is the mass of the DEF particles set with the parameter `Mass` in amu.


For this section there is no parameter order.

#### SU (to include a substrate)

This section aims at setting the substrate parameters.

- To set the thickness of the substrate you can use the `Thickness` parameter in Angstrom.
- The density is in particle/nm^3 and set using the `Density` parameter.
- The topology version is set using the `Version` parameter.
- The type of particle to use in the substrate is defined using the `SuType` parameter.

The following parameters are optional:

- In the case of membranes the bilayer height can be set relative to the substrate thickness using `BilayerHeight`.  in Angstrom.
- For rough or patterned substrates a parameter `RUGOSITY` can be set to describe the roughness. The `RUGOSITY` is described by a set of parameters defined with parameter:value and separated by | .

For this section there is no parameter order.

#### WALL (to include walls)

For the wall only the topology version is used as parameter with `Version`.

For this section there is no parameter order.

#### MONO (to add a monolayer)

This section aims at setting the monolayer parameters to define trilayers.

- The number of lipids is set using `NbLipidsM` or alternatively the area per lipid can be used to set the number of lipid using `APL` parameter (To implement)
- The height in the box of the monolayer is set using `LzM` in Angstrom.

For this section there is no parameter order.

#### PROTOCOL

This section sets the limit between system definition and simulation steps.

If walls, defo or substrate are set in the system, parameters `DEFO`, `WALL` and `SU` must be set in the same row as the `PROTOCOL` section.
These parameters define the protocol for the simulations steps and consists in parameter preset names (defined in GROMACS_Default/[SU,DEFO,WALL]) separated by +.
There should be as many name as simulations steps.

Example for 6 simulation steps:

<table style="width:100%">
  <tr>
    <th> PROTOCOL </th>
    <th> DEFO</th> 
    <td> d1+d4+d2+d3+d1+d6 </td>
    <th> WALL </th>
    <td> w2+w3+w1+w5+w1+w6 </td>
    <th> SU </th>
    <td> s1+s4+s2+s3+s1+s1 </td>
  </tr>
</table>


#### INIT or INPUT

The first step of the simulation protocol is the creation of the system. It can be either initialised using the `INIT` section (no parameters needed) or the `INPUT` that will use files from other simulations as starting point.

If `INPUT` is used the following parameters are needed:

- the path to a system configuration file is selected using the `GRO` parameter. If the .top, .tpr and .ndx have similar naming convention there is no need to specify them. Otherwise their path are set separately using the `TOP`, `TPR` and `NDX` parameters resp.
- If you want to add a substrate to the existing configuration you can select the substrate you want to use (path to a .gro file) using the `SU` parameter. The thickness, particle type and version of the substrate should be set in the SU section.
- If a DEFO is already in the input configuration, you need to specify the parameters of the DEFO. In particular, the version for the topology.
- If you want to generate a substrate having the same dimensions has the input, the parameter `GEN_SU` will allow you to do that. The thickness, particle type and version of the substrate to generate should be set in the SU section.
- Using the parameter `TRANSLATE` you can translate the system in each direction *X Y Z* in Angstrom. (gmx trjconv is used with pbc)
- With the parameter `NEWBOX` you can increase the box dimensions *LX LY LZ* in Angstrom.
- Using the `THRESHOLD-Z` along with `THRESHOLD_ATOMS` parameters you can remove particles of type `THRESHOLD_ATOMS` above a given z value in Angstrom set with `THRESHOLD-Z`.

#### Simulation steps

Each simulation step is described using two rows.

- The first row defines and the name of the step and the parameters to set.
- The second row sets the parameters values.

Example:
<table style="width:100%">
  <tr>
    <th>Simulation step 1 name </th>
    <th>parameter 1</th> 
    <th>Parameter 2 </th>
    <th>Parameter 3 </th>
    <th> ... </th>
  </tr>
  <tr>
    <td></td>
    <td>value 1</td> 
    <td>value 2</td>
    <td>value 3</td>
    <td> ... </td>
  </tr>
  <tr>
    <th>Simulation step 2 name </th>
    <th>parameter 1</th> 
    <th>Parameter 2 </th>
    <th>Parameter 3 </th>
    <th> ... </th>
  </tr>
  <tr>
    <td></td>
    <td>value 1</td> 
    <td>value 2</td>
    <td>value 3</td>
    <td> ... </td>
  </tr>
</table>
