from pathlib import Path
import csv
import json
import numpy as np
from Bio.PDB import MMCIFParser

RESULTS_DIR = Path(snakemake.input.boltz_results)
OUTPUT_CSV = Path(snakemake.output[0])


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def find_one(folder, pattern):
    files = sorted(folder.rglob(pattern))
    return files[0] if files else None


def get_manifest_id(manifest_path):
    return load_json(manifest_path)["records"][0]["id"]


def pair_iptm(confidence, scored_chain, aligned_chain):
    return confidence.get("pair_chains_iptm", {}).get(scored_chain, {}).get(aligned_chain, "")


def load_atoms_by_direct_chain(cif_path):
    structure = MMCIFParser(QUIET=True).get_structure("model", cif_path)
    atoms_by_chain = {}

    for atom in structure.get_atoms():
        chain_id = atom.get_parent().get_parent().id
        atoms_by_chain.setdefault(chain_id, []).append(atom)

    return atoms_by_chain


def count_contacts(candidate_atoms, ligand_atoms, cutoff=4.0):
    if not candidate_atoms or not ligand_atoms:
        return 0

    c = np.array([x.coord for x in candidate_atoms])
    l = np.array([x.coord for x in ligand_atoms])
    
    # matrice des distances de paires vectorisée
    d = np.linalg.norm(c[:, None, :] - l[None, :, :], axis=2)
    return int(np.sum(d <= cutoff))


def count_interface_residues(candidate_atoms, ligand_atoms, cutoff=4.0):
    if not candidate_atoms or not ligand_atoms:
        return 0

    c = np.array([x.coord for x in candidate_atoms])
    l = np.array([x.coord for x in ligand_atoms])
    
    d = np.linalg.norm(c[:, None, :] - l[None, :, :], axis=2)
    candidate_indices_under_cutoff = np.any(d <= cutoff, axis=1)
    
    residues = set()
    for idx, under_cutoff in enumerate(candidate_indices_under_cutoff):
        if under_cutoff:
            atom = candidate_atoms[idx]
            residue = atom.get_parent()
            residues.add((residue.get_parent().id, residue.id[1], residue.resname))

    return len(residues)


def min_distance(atoms_a, atoms_b):
    if not atoms_a or not atoms_b:
        return ""

    a = np.array([x.coord for x in atoms_a])
    b = np.array([x.coord for x in atoms_b])
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)

    return float(np.min(d))


def ligand_coordination(candidate_atoms, ligand_atoms, cutoff=3.0):
    if not candidate_atoms or not ligand_atoms:
        return 0, ""

    ligand_coords = np.array([a.coord for a in ligand_atoms])
    residue_min_distances = {}

    for atom in candidate_atoms:
        # seuls les atomes possédant des doublets électroniques libres..
        if atom.element not in {"N", "O", "S", "Se"}:
            continue

        d = float(np.min(np.linalg.norm(ligand_coords - atom.coord, axis=1)))

        if d <= cutoff:
            residue = atom.get_parent()
            res_id = (residue.get_parent().id, residue.id[1])
            
            if res_id not in residue_min_distances or d < residue_min_distances[res_id]:
                residue_min_distances[res_id] = d

    if not residue_min_distances:
        return 0, ""

    coordinating_distances = list(residue_min_distances.values())
    return (
        len(coordinating_distances),
        float(np.mean(coordinating_distances)),
    )


def analyze_candidate(candidate_dir):
    manifest = find_one(candidate_dir, "manifest.json")
    confidence_json = find_one(candidate_dir, "confidence_*_model_0.json")
    cif = find_one(candidate_dir, "*_model_0.cif")

    if not manifest or not confidence_json or not cif:
        print(f"Skip {candidate_dir.name}: fichiers Boltz manquants")
        return None

    candidate_id = get_manifest_id(manifest)
    confidence = load_json(confidence_json)

    atoms = load_atoms_by_direct_chain(cif)
    candidate_atoms = atoms.get("A", [])
    l1_atoms = atoms.get("L1", [])
    l2_atoms = atoms.get("L2", [])
    
    if not candidate_atoms:
        print(f"Warning {candidate_dir.name}: Chaine de protéine 'A' vide.")

    l2_coord_count, l2_coord_mean_dist = ligand_coordination(candidate_atoms, l2_atoms)

    #l1_iptm = pair_iptm(confidence, "1", "0") or pair_iptm(confidence, "L1", "A")
    #l2_iptm = pair_iptm(confidence, "2", "0") or pair_iptm(confidence, "L2", "A")
    record = load_json(manifest)["records"][0]
    chain_ids = {
        chain["chain_name"]: str(chain["chain_id"])
        for chain in record["chains"]
    }

    l1_iptm = pair_iptm(confidence, chain_ids["L1"], chain_ids["A"])
    l2_iptm = pair_iptm(confidence, chain_ids["L2"], chain_ids["A"])

    return {
        "candidate_id": candidate_id,
        "confidence_score": confidence.get("confidence_score", ""),
        "ligand_iptm": confidence.get("ligand_iptm", ""),
        "complex_iplddt" : confidence.get("complex_iplddt", ""),
        "L1_candidate_iptm": l1_iptm,
        "L2_candidate_iptm": l2_iptm,
        "L1_candidate_contacts_4A": count_contacts(candidate_atoms, l1_atoms),
        "L2_candidate_contacts_4A": count_contacts(candidate_atoms, l2_atoms),
        "L1_candidate_interface_residues": count_interface_residues(candidate_atoms, l1_atoms),
        "L2_candidate_interface_residues": count_interface_residues(candidate_atoms, l2_atoms),
        "L1_L2_min_distance": min_distance(l1_atoms, l2_atoms),
        "L2_coordination_count_3A": l2_coord_count,
        "L2_mean_coordination_distance_3A": l2_coord_mean_dist,
    }


rows = []
for candidate_dir in sorted(p for p in RESULTS_DIR.iterdir() if p.is_dir()):
    row = analyze_candidate(candidate_dir)
    if row:
        rows.append(row)

fieldnames = [
    "candidate_id", "confidence_score", "ligand_iptm", "complex_iplddt", "L1_candidate_iptm", "L2_candidate_iptm",
    "L1_candidate_contacts_4A", "L2_candidate_contacts_4A", "L1_candidate_interface_residues",
    "L2_candidate_interface_residues", "L1_L2_min_distance", "L2_coordination_count_3A",
    "L2_mean_coordination_distance_3A"
]

OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
