'''
Helper script to convert Abver incprd coordinate format to start/finish format.
'''
import sys

def split_coords(infile, outfile):
    with open(infile, 'r') as inf, open(outfile, 'w') as outf:
        lines = inf.readlines()

        if len(lines) < 2:
            raise ValueError("Input file must have at least two lines (header + atom count).")

        # ignore first line, read second for atom count (not strictly required but validated)
        atom_count_line = lines[1].strip()
        try:
            atom_count = int(atom_count_line)
        except ValueError:
            raise ValueError("Second line must contain number of atoms as an integer.")

        # process remaining lines (from third line onward)
        data_lines = lines[2:]
        coords = []
        for i, line in enumerate(data_lines, start=3):
            parts = line.split()
            if not parts:
                continue
            if len(parts) % 3 != 0:
                raise ValueError(f"Line {i} does not contain coordinates in multiples of 3: {line!r}")
            # split into groups of 3 (x,y,z)
            for j in range(0, len(parts), 3):
                x, y, z = parts[j:j+3]
                coords.append((x, y, z))

        if len(coords) != atom_count:
            # warn but still write coordinates; you can change to raise if strict enforcement desired
            print(f"Warning: parsed {len(coords)} coordinates but second line reported {atom_count} atoms.", file=sys.stderr)

        # write each atom on its own line as "x y z"
        for x, y, z in coords:
            outf.write(f"{x}\t{y}\t{z}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python split_coords.py coords.inpcrd start")
        sys.exit(1)
    split_coords(sys.argv[1], sys.argv[2])
