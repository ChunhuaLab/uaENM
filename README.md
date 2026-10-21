# uaENM: Integrating United-atom Modeling with a Corresponding Potential Elastic Network Model for Characterization of Protein Dynamics

uaENM combines united-atom protein modeling with an elastic network potential to calculate residue fluctuations.

Authors: Xinyu Zhang, Zhixiang Wu, Jilong Zhang, Jingjie Su, Long Zhao, Chunhua Li.

## Table of Contents

- [Requirements](#requirements)
- [Project Structure](#project-structure)
- [Full Pipeline](#full-pipeline)
- [Calculation Method](#calculation-method)
- [Help](#help)

## Requirements

- MATLAB R2020a or later. No additional MATLAB toolbox is required.
- The SPICA model conversion tool, available at [SPICA](https://www.spica-ff.org), for preparing the united-atom structure and topology.

## Project Structure

```text
uaENM_code/
|-- README.md
|-- codes/
|   |-- res_type_spica.m         # SPICA force-field parameters
|   |-- read_top.m               # Structure and topology reader
|   |-- pdbread_uaENM.m          # United-atom PDB parser
|   |-- hessian_2.m              # Two-body Hessian
|   |-- hessian_3.m              # Bond-angle Hessian
|   |-- hessian_4.m              # Dihedral Hessian
|   |-- confirm_twobody_pir.m    # Non-bonded pair selection
|   |-- calculate_uaENM.m        # Residue fluctuations and B-factor fitting
|   `-- main.m                   # Run parameter optimization and save results
`-- example/
    |-- 3mp9-y.pdb
    |-- 3mp9-cg.pdb
    `-- 3mp9.top
```

## Full Pipeline

### 1. Prepare the united-atom structure and topology

The complete 3MP9 example is included. For another protein, replace the three input files in `example/` with matching `<id>-y.pdb`, `<id>-cg.pdb` and `<id>.top` files. Keep one protein in this folder per run.

Prepare the input files with the SPICA conversion tool:

```bash
cg_spica map2cg 3mp9-y.pdb 3mp9-cg.pdb
cg_spica ENM 3mp9-cg.pdb 3mp9.top -aapdb 3mp9-y.pdb -dssp /usr/bin/dssp -pspica
```

The reference `<id>-y.pdb` must contain experimental B-factors in columns 61-66 of its ATOM CA records. Its residues must correspond to the CG backbone sites in the same order. Use a single PDB model without duplicate CA alternate locations.

### 2. Run calculation

Open `codes/main.m` in MATLAB and click **Run**. The script searches all parameter combinations, selects the best eligible Pearson PCC, fits the prediction amplitude to the experimental B-factors and displays the comparison figure. No separate launcher or arguments are needed.

All results are saved directly in `example/`:

| File | Content |
|---|---|
| `resul.txt` | Best PCC, cutoff and nrc (`n_x`) |
| `uaENM_result.png` | Experimental and fitted uaENM B-factor curves |
| `msf.txt` | Raw mean-square fluctuation for each residue, one value per line |
| `Bfactor.txt` | Fitted B-factor for each residue, one value per line |

The best PCC and parameters also remain in the MATLAB workspace as `Smax_uaENM`, `cutoff` and `R_ab`.

## Calculation Method

- Search grid: cutoff `10:20` Å and nrc `1.5:0.5:4.0`, giving 66 combinations.
- Non-bonded interactions use the current vdW potential without electrostatics. Spring constants are capped at 10; bonded, bond-angle and dihedral terms retain their current definitions.
- Eigenvalues and eigenvectors are sorted in ascending order. Non-mass-weighted CA fluctuations use modes 7 through the final mode, with no extra eigenvalue filtering or normalization:

  `MSF_i = sum_j ((v_ix,j^2 + v_iy,j^2 + v_iz,j^2) / lambda_j), j = 7,...,3N`.

- Eligible combinations require `lambda7 >= 1e-4` and a raw capped non-bonded spring fraction `< 0.05`. The maximum **unrounded** Pearson PCC is selected. Exact ties use the smaller cutoff, then the smaller nrc.
- The comparison curve uses an ordinary least-squares fit with an intercept, `B_fitted = a + b * MSF`, to adjust both amplitude and offset. The plot compares two B-factor curves in Å²; `msf.txt` retains the raw model output without an additional temperature factor.


## Help

For questions, contact chunhuali@bjut.edu.cn.
