import itertools
import os
import tempfile

import click
import numpy as np
import pandas as pd
from Bio.PDB import PDBList, PDBParser


def get_residue_id(residue):
    _, id, _ = residue.get_id()
    return id


@click.command()
@click.option("--pdb-code")
@click.option("--pdb-file")
@click.option("--cutoff", default=6)
def cli(pdb_code, pdb_file, cutoff):
    # Get structure
    if pdb_code is not None:
        pdb = PDBList()
        tmp_dir = tempfile.mkdtemp()
        pdb_file = pdb.retrieve_pdb_file(pdb_code, pdir=tmp_dir, file_format="pdb")
        os.rmdir("obsolete")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # Calculate the number of pairs
    residue_list = [
        r
        for r in structure.get_residues()
        if r.get_resname() not in ["WAT", "Na+", "Cl-", "HOH"]
    ]
    pair_list = list(itertools.product(residue_list, residue_list))
    sides = ["left", "right"]
    pair_df = pd.DataFrame(data=pair_list, columns=sides)
    for side in sides:
        pair_df[f"name_{side}"] = pair_df[side].apply(lambda r: r.get_resname())
        pair_df[f"id_{side}"] = pair_df[side].apply(get_residue_id)
        pair_df[f"center_of_mass_{side}"] = pair_df[side].apply(
            lambda r: r.center_of_mass()
        )
        pair_df[f"center_of_geometry_{side}"] = pair_df[side].apply(
            lambda r: r.center_of_mass(geometric=True)
        )
    # "center" is usually called "center of geometry"
    # Define atom coordinates as A_i, the atoms belong to a residue that contain n atoms
    # The center of geometry of a residue, R is:
    # R = \frac{\Sigma_i A_i}{n}
    for method in ["center_of_geometry", "center_of_mass"]:
        diff = pair_df[f"{method}_left"] - pair_df[f"{method}_right"]
        pair_df[f"distance_{method}"] = (diff * diff).apply(
            lambda diff2: np.sqrt(sum(diff2))
        )

    pair_df["distance_nearest"] = pair_df.apply(
        lambda row: min(
            abs(
                np.diff(
                    list(
                        itertools.product(
                            row["left"].get_atoms(), row["right"].get_atoms()
                        )
                    )
                )
            )
        )[0],
        axis=1,
    )

    for side in sides:
        pair_df = pair_df.drop(
            columns=[side, f"center_of_mass_{side}", f"center_of_geometry_{side}"]
        )

    pair_df.to_csv(f"{os.path.basename(pdb_file)}.csv", index=False)

    npairs_nearest = pair_df[pair_df.distance_nearest < cutoff].distance_nearest.count()
    npairs_CoM = pair_df[
        pair_df.distance_center_of_mass < cutoff
    ].distance_center_of_mass.count()
    npairs_CoG = pair_df[
        pair_df.distance_center_of_geometry < cutoff
    ].distance_center_of_geometry.count()
    print(
        f"""
Cutoff Length: {cutoff} Angstrom
Distances between:
\tNearest: {npairs_nearest}
\tCenter of mass: {npairs_CoM}
\tCenter of geometry: {npairs_CoG}
"""
    )


if __name__ == "__main__":
    cli()
