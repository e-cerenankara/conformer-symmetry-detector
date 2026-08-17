import argparse
import math


def parse_xyz_frames(file):
    # reads an xyz file, a file can have one or more frames in it
    with open(file) as f:
        lines = f.readlines()

    frames = []
    index = 0
    while index < len(lines):
        n_atoms = int(lines[index])
        atoms = []
        for atom_index in range(n_atoms):
            parts = lines[index + 2 + atom_index].split()
            symbol = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            atoms.append((symbol, x, y, z))
        frames.append(atoms)
        index += n_atoms + 2

    return frames


def load_frames(files):
    # loads frames from one or more files and labels each one
    labeled_frames = []
    for file in files:
        frames = parse_xyz_frames(file)
        for i, atoms in enumerate(frames):
            if len(frames) == 1:
                label = file
            else:
                label = f"{file} (frame {i + 1})"
            labeled_frames.append((label, atoms))
    return labeled_frames


def same_atom_types(atoms1, atoms2):
    # checks if both frames have the same elements
    types1 = sorted(atom[0] for atom in atoms1)
    types2 = sorted(atom[0] for atom in atoms2)
    return types1 == types2


def pairwise_distances(atoms):
    # all pairwise distances between atoms, sorted
    distances = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            _, x1, y1, z1 = atoms[i]
            _, x2, y2, z2 = atoms[j]
            dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
            distances.append(math.sqrt(dx * dx + dy * dy + dz * dz))
    distances.sort()
    return distances


def same_distances(distances1, distances2, tolerance=1e-4):
    # compares with a small tolerance instead of exact equality
    if len(distances1) != len(distances2):
        return False
    for d1, d2 in zip(distances1, distances2):
        if abs(d1 - d2) > tolerance:
            return False
    return True


def group_identical_frames(labeled_frames):
    # groups frames that turn out to be the same conformation
    groups = []
    used = set()

    for i in range(len(labeled_frames)):
        if i in used:
            continue

        label_i, atoms_i = labeled_frames[i]
        group = [label_i]
        distances_i = pairwise_distances(atoms_i)

        for j in range(i + 1, len(labeled_frames)):
            if j in used:
                continue
            label_j, atoms_j = labeled_frames[j]
            if same_atom_types(atoms_i, atoms_j):
                distances_j = pairwise_distances(atoms_j)
                if same_distances(distances_i, distances_j):
                    group.append(label_j)
                    used.add(j)

        groups.append(group)

    return groups


def print_report(labeled_frames, groups):
    print("Symmetry Detection Results:")
    print("-" * 50)
    print("Total frames:", len(labeled_frames))
    print("Unique conformations:", len(groups))
    print()

    for group in groups:
        if len(group) > 1:
            print("These frames are identical:")
            for label in group:
                print(" -", label)
        else:
            print("Unique:", group[0])
        print()
    print("-" * 50)


def main():
    parser = argparse.ArgumentParser(description="Detect identical conformations across one or more XYZ files.")
    parser.add_argument("files", nargs="+", help="Path(s) to .xyz file(s)")
    args = parser.parse_args()

    labeled_frames = load_frames(args.files)
    groups = group_identical_frames(labeled_frames)
    print_report(labeled_frames, groups)


if __name__ == "__main__":
    main()