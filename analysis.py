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
    """Return (center [x,y,z], size [x,y,z]) from data/box.yaml, using CONFIG for overrides."""
    box_path = config.get("box")
    if not box_path:
        raise ValueError("config must provide a 'box' path")

    with open(box_path, encoding="utf-8") as box_file:
        box_data = yaml.safe_load(box_file) or {}
    if not isinstance(box_data, dict):
        raise ValueError(f"Box file {box_path!r} must contain a YAML mapping")

    center = config.get("center")
    if center is None:
        center = box_data.get("center")
    size = config.get("size")
    if size is None:
        size = box_data.get("size")
    if center is None or size is None:
        raise ValueError(
            f"Box center and size are missing; prepare receptor data and populate {box_path!r}"
        )

    try:
        center = [float(value) for value in center]
        size = [float(value) for value in size]
    except (TypeError, ValueError) as error:
        raise ValueError("Box center and size must each contain three numbers") from error
    if (len(center) != 3 or len(size) != 3
            or not np.isfinite(center).all() or not np.isfinite(size).all()
            or any(value <= 0 for value in size)):
        raise ValueError(
            "Box center and size must each contain three finite values; sizes must be positive"
        )

    return center, size


# ---------- Student A ----------
def dock_ligand(receptor, ligand, center, size, exhaustiveness):
    """Dock one ligand with Vina. Return (best_score, pose_pdbqt_string).
    Hint: v = Vina(sf_name='vina'); v.set_receptor(...); v.set_ligand_from_file(...);
    v.compute_vina_maps(center=..., box_size=...); v.dock(exhaustiveness=..., n_poses=5)"""
    raise NotImplementedError


def pose_rmsd(pose_pdbqt, crystal_pdbqt):
    """Heavy-atom RMSD between the docked pose and the crystal ligand (same atom order
    is NOT guaranteed - match by element and nearest neighbour, or use the symmetry-
    corrected RMSD from meeko/rdkit)."""
    raise NotImplementedError


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


def main():
    center, size = load_box(CONFIG)
    exhaustiveness = int(CONFIG.get("exhaustiveness", 8))
    if exhaustiveness < 1:
        raise ValueError("Set 'exhaustiveness' to a positive integer in config.yaml")

    receptor = CONFIG["receptor"]
    ligand_dir = CONFIG["ligand_dir"]
    crystal_ligand = CONFIG["crystal_ligand"]
    sotorasib = os.path.join(ligand_dir, "sotorasib.pdbqt")

    redock_score, pose_pdbqt = dock_ligand(receptor, sotorasib, center, size, exhaustiveness)
    with open(crystal_ligand, encoding="utf-8") as crystal_file:
        crystal_pdbqt = crystal_file.read()
    rmsd = pose_rmsd(pose_pdbqt, crystal_pdbqt)

    os.makedirs("results", exist_ok=True)
    pd.DataFrame([{"ligand": "sotorasib", "score": redock_score, "rmsd": rmsd}]).to_csv(
        "results/redock.csv", index=False
    )

    scores = dock_all(
        receptor, ligand_dir, center, size, exhaustiveness,
        exclude_ligands={"sotorasib"},
    )
    scores = pd.concat(
        [scores, pd.DataFrame([{"ligand": "sotorasib", "score": redock_score}])],
        ignore_index=True,
    )
    plot_ranking(scores)
    print(summary_sentence(rmsd, scores))


if __name__ == "__main__":
    main()
