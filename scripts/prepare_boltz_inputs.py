from pathlib import Path
import re
import yaml

hits_fasta = snakemake.input.hits
ligands_yaml = snakemake.input.ligands
output_dir = Path(snakemake.output[0])

output_dir.mkdir(parents=True, exist_ok=True)


def parse_fasta(fasta_path):
    records = []
    current_id = None
    current_seq = []

    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_seq)))

                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        records.append((current_id, "".join(current_seq)))

    return records


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


with open(ligands_yaml, "r", encoding="utf-8") as f:
    ligand_data = yaml.safe_load(f)

ligands = ligand_data.get("ligands", [])
if not ligands:
    raise ValueError("aucun ligand trouvé dans le fichier origin_ligands.yaml")

hits = parse_fasta(hits_fasta)
if not hits:
    raise ValueError("aucune séquence trouvée dans le fichier FASTA des hits")

for hit_id, sequence in hits:
    ligand_blocks = []

    for i, ligand in enumerate(ligands, start=1):
        smiles = ligand.get("smiles")

        if not smiles:
            raise ValueError(
                f"smiles manquant pour {ligand.get('name')} ({ligand.get('chebi_id')})"
            )

        ligand_blocks.append(
            {
                "ligand": {
                    "id": f"L{i}",
                    "smiles": smiles,
                }
            }
        )

    job_name = safe_name(hit_id)

    boltz_input = {
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": sequence,
                }
            },
            *ligand_blocks,
        ],
    }

    output_file = output_dir / f"{job_name}.yaml"

    with open(output_file, "w", encoding="utf-8") as out:
        yaml.safe_dump(
            boltz_input,
            out,
            sort_keys=False,
            default_flow_style=False,
        )
