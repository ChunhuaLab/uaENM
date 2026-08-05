# uaENM:Integrating United-atom Modeling with Corresponding Potential Elastic Network Model for Characterization of Protein Dynamics

We propose a united atom ENM (uaENM), which incorporates a united-atom protein modeling and an enhanced potential for for high-throughput Characterization of protein dynamic information.

Authors: Xinyu Zhang, Zhixiang Wu, Jilong Zhang, Jingjie Su, Long Zhao, Chunhua Li.

---

## Table of Contents

- [Requirements](#requirements)
- [Project Structure](#project-structure)
- [Full Pipeline](#full-pipeline)
- [Step 1: Prepare the united atom structure and topology for a protein](#step-1-Prepare-the-united-atom-structure-and-topology-for-a-protein)
- [Step 2: Run calculation](#step-2-Run-calculation)

---

## Requirements

> The protein model conversion tool provided by the SPICA model (https://www.spica-ff.org) for generating the united atom topology for protein (For example file 3pe9-cg.pdb and 3pe9.top). 

>MATLAB R2020+, Parallel Computing Toolbox.

---

## Project Structure

```
uaENM/
©À©¤©¤ README.md
©À©¤©¤ codes/
©¦   ©À©¤©¤ res_type_spica.m         # SPICA Force Field Parameter Repository 
©¦   ©À©¤©¤ read_top.m               # United-Atom Protein Model Structure Organizer
©¦   ©À©¤©¤ pdbread_uaENM.m          # United-Atom Protein Model Structure Parser
©¦   ©À©¤©¤ hessian_2.m              # Two-Body Potential Hessian Matrix Calculator
©¦   ©À©¤©¤ hessian_3.m              # Three-Body Potential Hessian Matrix Calculator
©¦   ©À©¤©¤ hessian_4.m              # Four-Body Potential Hessian Matrix Calculator
©¦   ©À©¤©¤ confirm_twobody_pir.m    # Non-Bonded Two-Body Atomic Pair Determination
©¦   ©À©¤©¤ calculate_uaENM.m        # uaENM calculation Main Program
©¦   ©¸©¤©¤ main.m                   # Parameter Optimization and Result Analysis
©¸©¤©¤ example/
    ©À©¤©¤ 3pe9-y.pdb
    ©À©¤©¤ 3pe9-cg.pdb
    ©À©¤©¤ 3pe9.top
```

## Full Pipeline

### 1. Prepare the united atom structure and topology for a protein.

Note : First, get the conversion tool on https://www.spica-ff.org
An example of a generated structure:
```
cg_spica map2cg 3pe9-y.pdb 3pe9-cg.pdb

```

An example of a generated topology :

```

 cg_spica ENM 3pe9-cg.pdb 3pe9.top -aapdb 3pe9-y.pdb -dssp /usr/bin/dssp -pspica

```


### 2. Run calculation

After running the main.m program, the fluctuation of the protein residue, the optimal parameters and the comparison image with the experimental B-factor will be obtained.

## Help

For any questions, please contact us by chunhuali@bjut.edu.cn.
