# Save the trajectory frame closest to the average structure

This VMD Tcl script loads a selected set of Amber NetCDF trajectories, aligns every loaded frame to frame 0 using the requested protein backbone residues, calculates the average coordinates for those residues, and saves the loaded frame with the lowest RMSD to that average as a PDB.

The saved PDB is an actual trajectory frame. It is not a structure made from averaged coordinates.

## Requirements

- VMD with support for Amber `parm7` topology files and NetCDF trajectories
- One topology file and all requested trajectory files in the directory from which VMD is launched
- Matching atom order and compatible topology across all trajectory files

## Input filenames

The script expects these names:

```text
<system>_md.prmtop
<system>_step<step><replicate>.nc
```

For example, system `beta2_tyrosine`, steps `7-10`, and replicates `a-d` require `beta2_tyrosine_md.prmtop` and trajectory files such as `beta2_tyrosine_step7a.nc` through `beta2_tyrosine_step10d.nc`.

All requested step/replicate combinations must exist. The script loads them in ascending step order, then ascending replicate order. It starts at frame 0 in each file and uses the stride entered at the prompt.

## Run

Open a terminal in the directory containing the topology and trajectory files, then run:

```bash
vmd -dispdev text -e /path/to/save_average_frame.tcl
```

Enter the requested values when prompted:

1. **System name:** the shared filename prefix, such as `beta2_tyrosine`
2. **Steps:** one step or an inclusive range, such as `7-10`
3. **Replicates:** one lowercase letter or an inclusive letter range, such as `a-d`
4. **Stride:** a positive integer, such as `100`
5. **Alignment residue ranges:** residue numbers and ranges separated by spaces

For the last prompt, choose the most stable regions of the protein for alignment. For a GPCR, the transmembrane domains are a common choice; avoid flexible loops or tails where possible. Enter residue numbers only, without the word `resid` or parentheses. For example:

```text
37 to 59 66 to 90 108 to 137 146 to 164 203 to 230 266 to 291 311 to 326
```

The script converts that to this VMD selection:

```text
backbone and (resid 37 to 59 66 to 90 108 to 137 146 to 164 203 to 230 266 to 291 311 to 326)
```

Use residue numbering as it appears in the loaded topology. The same selection is used for alignment and RMSD comparison.

## Output

The script writes `<system>_avg_frame.pdb` in the current working directory. If that file already exists, it asks before replacing it. The reported frame index is VMD's zero-based index in the combined, strided trajectory.

The included `.gitignore` excludes generated representative-frame PDBs and the system-specific topology and trajectory files, which are typically large and should not be committed with the script.

## Method

Each loaded frame is fitted to frame 0 using the selected backbone atoms, and the same transformation is applied to all atoms in that frame. The script averages the selected atom coordinates over the aligned frames, compares each aligned frame with that average by RMSD, and writes all atoms from the closest frame.

Because an average structure can fall between conformations sampled by the trajectory, the closest loaded frame is a representative snapshot rather than the coordinate average itself.
