# Conformer Symmetry Detection

A Python script that checks if different XYZ structures are actually the same conformation, even when moved or rotated in space.

## What is a conformer?

A conformer is a shape a molecule can take by rotating around a single bond. Ethane has two well known ones. One of them is staggered, where the hydrogens stay as far apart as possible, and eclipsed, where they line up directly across from each other.

## How this works

Instead of comparing raw coordinates, the script compares the distances between all pairs of atoms. Two copies of the same conformer can have totally different coordinates if one is moved or rotated, but the distances between their atoms stay the same. Comparing distances instead of coordinates makes the check work no matter how the structure is placed.

1. Read each XYZ file.
2. Check that two frames have the same elements before comparing anything else.
3. Calculate every pairwise distance in a frame, and sort them.
4. Compare two frames by checking if their distance lists match within a small tolerance.
5. Group frames that match into the same conformation.

## Requirements

Python 3, standard library only.

## Installation

```bash
git clone https://github.com/e-cerenankara/conformer-symmetry.git
cd conformer-symmetry
```

## Usage

```bash
python symmetry_detector.py file1.xyz file2.xyz file3.xyz
```

## Example

```bash
python symmetry_detector.py examples/ethane_a.xyz examples/ethane_b.xyz examples/ethane_c.xyz
```

Output

```
Symmetry Detection Results:
--------------------------------------------------
Total frames: 3
Unique conformations: 2

These frames are identical:
 - examples/ethane_a.xyz
 - examples/ethane_b.xyz

Unique: examples/ethane_c.xyz
--------------------------------------------------
```

`ethane_a.xyz` and `ethane_b.xyz` are the same staggered conformer, just placed differently in space. `ethane_c.xyz` is eclipsed, a genuinely different shape, so it stays on its own.


