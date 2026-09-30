"""
Shared analysis script. Both students edit THIS file on their own branches.

Student A implements : dock_ligand, pose_rmsd
Student B implements : dock_all, plot_ranking
BOTH implement       : load_box, main, summary_sentence   <- expect a merge conflict here

Run: python analysis.py
"""
import yaml
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONFIG = yaml.safe_load(open("config.yaml"))

import glob, os
import numpy as np
from vina import Vina

# ---------- BOTH ----------
def load_box(config):
    """Return docking box center and size."""

    with open(config["box"]) as f:
        box = yaml.safe_load(f)

    center = config.get("center", box["center"])
    size = config.get("size", box["size"])

    return center, size


# ---------- Student A ----------
def dock_ligand(receptor, ligand, center, size, exhaustiveness):
    """Dock one ligand with Vina. Return (best_score, pose_pdbqt_string)."""
    v = Vina(sf_name="vina")

    v.set_receptor(receptor)
    v.set_ligand_from_file(ligand)

    v.compute_vina_maps(
        center=center,
        box_size=size
    )

    v.dock(
        exhaustiveness=exhaustiveness,
        n_poses=5
    )

    best_score = float(v.energies(n_poses=1)[0][0])
    pose = v.poses(n_poses=1)

    return best_score, pose


def pose_rmsd(pose_pdbqt, crystal_pdbqt):
    """Calculate heavy-atom RMSD between two PDBQT poses."""

    def read_atoms(pdbqt):
        atoms = []

        for line in pdbqt.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                atom_type = line[77:79].strip().upper()

                if atom_type.startswith("H"):
                    continue

                element = atom_type[0]

                atoms.append({
                    "element": element,
                    "coord": np.array([x, y, z], dtype=float)
                })

        return atoms

    pose_atoms = read_atoms(pose_pdbqt)
    crystal_atoms = read_atoms(crystal_pdbqt)

    if len(pose_atoms) != len(crystal_atoms):
        raise ValueError(
            "Docked and crystal poses have different numbers of heavy atoms."
        )

    matched_pose = []
    matched_crystal = []
    used = set()

    for crystal_atom in crystal_atoms:
        candidates = [
            (i, atom)
            for i, atom in enumerate(pose_atoms)
            if i not in used
            and atom["element"] == crystal_atom["element"]
        ]

        if not candidates:
            raise ValueError(
                f"No matching atom found for element {crystal_atom['element']}."
            )

        best_index, best_atom = min(
            candidates,
            key=lambda item: np.linalg.norm(
                item[1]["coord"] - crystal_atom["coord"]
            )
        )

        used.add(best_index)
        matched_crystal.append(crystal_atom["coord"])
        matched_pose.append(best_atom["coord"])

    P = np.array(matched_pose)
    Q = np.array(matched_crystal)

    P_centered = P - P.mean(axis=0)
    Q_centered = Q - Q.mean(axis=0)

    H = P_centered.T @ Q_centered
    U, _, Vt = np.linalg.svd(H)

    R = Vt.T @ U.T

    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    P_aligned = P_centered @ R

    rmsd = np.sqrt(
        np.mean(np.sum((P_aligned - Q_centered) ** 2, axis=1))
    )

    return float(rmsd)


# ---------- Student B ----------
def dock_all(receptor, ligand_dir, center, size, exhaustiveness, exclude_ligands=()):
    """Dock ligand files, optionally excluding basenames; return ligand and score columns."""
    ligand_files = sorted(glob.glob(os.path.join(ligand_dir, "*.pdbqt")))
    if not ligand_files:
        raise FileNotFoundError(f"No .pdbqt ligands found in {ligand_dir!r}")
    excluded = {str(name).lower() for name in exclude_ligands}
    ligand_files = [
        path for path in ligand_files
        if os.path.splitext(os.path.basename(path))[0].lower() not in excluded
    ]
    if not ligand_files:
        raise ValueError("No ligands remain after applying exclude_ligands")

    vina = Vina(sf_name="vina")
    vina.set_receptor(receptor)
    vina.compute_vina_maps(center=center, box_size=size)

    results = []
    for ligand_file in ligand_files:
        vina.set_ligand_from_file(ligand_file)
        vina.dock(exhaustiveness=exhaustiveness, n_poses=5)

        best_score = float(vina.energies(n_poses=1)[0][0])
        results.append({"ligand": os.path.splitext(os.path.basename(ligand_file))[0],
                        "score": best_score})

    return pd.DataFrame(results, columns=["ligand", "score"])


def plot_ranking(scores, out="results/ranking.png"):
    """Bar plot of Vina scores sorted best-first, sotorasib highlighted."""
    if not {"ligand", "score"}.issubset(scores.columns):
        raise ValueError("scores must contain 'ligand' and 'score' columns")
    if scores.empty:
        raise ValueError("scores must contain at least one ligand")

    ranking = scores.sort_values("score", ascending=True)
    ranking["score"] = pd.to_numeric(ranking["score"], errors="raise")
    if not np.isfinite(ranking["score"]).all():
        raise ValueError("scores must contain only finite numeric values")
    colors = ["#d1495b" if "sotorasib" in str(ligand).lower() else "#4c78a8"
              for ligand in ranking["ligand"]]

    figure, axis = plt.subplots(figsize=(max(7, len(ranking) * 1.1), 4.5))
    axis.bar(ranking["ligand"], ranking["score"], color=colors)
    axis.set_xlabel("Ligand")
    axis.set_ylabel("Vina score (kcal/mol)")
    axis.set_title("Docking score ranking (best first)")
    axis.tick_params(axis="x", labelrotation=35)
    figure.tight_layout()

    output_dir = os.path.dirname(os.path.abspath(out))
    os.makedirs(output_dir, exist_ok=True)
    figure.savefig(out, dpi=150)
    plt.close(figure)


# ---------- BOTH ----------
def summary_sentence(rmsd, scores):
    """One sentence: redocking RMSD, and the rank of sotorasib among the six ligands."""
    if not {"ligand", "score"}.issubset(scores.columns):
        raise ValueError("scores must contain 'ligand' and 'score' columns")

    sotorasib_mask = scores["ligand"].astype(str).str.contains("sotorasib", case=False)
    if not sotorasib_mask.any():
        raise ValueError("scores do not contain a sotorasib ligand")

    numeric_scores = pd.to_numeric(scores["score"], errors="raise")
    if not np.isfinite(numeric_scores).all():
        raise ValueError("scores must contain only finite numeric values")
    numeric_rmsd = float(rmsd)
    if not np.isfinite(numeric_rmsd) or numeric_rmsd < 0:
        raise ValueError("rmsd must be a finite, non-negative number")

    best_sotorasib_score = numeric_scores[sotorasib_mask].min()
    rank = int((numeric_scores < best_sotorasib_score).sum()) + 1
    return (
        f"The redocking RMSD was {numeric_rmsd:.2f} A, and sotorasib ranked "
        f"{rank} of {len(scores)} ligands by Vina score."
    )

    if sotorasib_rows.empty:
        raise ValueError("Could not find sotorasib in the docking results.")

    rank = int(sotorasib_rows.index[0]) + 1

    return (
        f"Sotorasib redocking RMSD was {rmsd:.2f} Å, "
        f"and it ranked {rank} out of {len(ordered)} ligands by Vina score."
    )
def main():
    center, size = load_box(CONFIG)

    receptor = CONFIG["receptor"]
    crystal_ligand = CONFIG["crystal_ligand"]
    exhaustiveness = CONFIG["exhaustiveness"]

    ligand_files = glob.glob(
        os.path.join(CONFIG["ligand_dir"], "*.pdbqt")
    )

    sotorasib = None

    for ligand in ligand_files:
        if "sotorasib" in os.path.basename(ligand).lower():
            sotorasib = ligand
            break

    if sotorasib is None:
        raise FileNotFoundError(
            "Could not find the sotorasib ligand in data/ligands."
        )

    score, pose = dock_ligand(
        receptor,
        sotorasib,
        center,
        size,
        exhaustiveness
    )

    rmsd = pose_rmsd(
        pose,
        open(crystal_ligand).read()
    )

    os.makedirs("results", exist_ok=True)

    pd.DataFrame([
        {
            "ligand": "sotorasib",
            "score": score,
            "rmsd": rmsd
        }
    ]).to_csv(
        "results/redock.csv",
        index=False
    )

    print(f"Sotorasib docking score: {score:.2f} kcal/mol")
    print(f"Sotorasib redocking RMSD: {rmsd:.2f} Å")


if __name__ == "__main__":
    main()