import re
import requests
import yaml

fasta_path = snakemake.input[0]
output_yaml = snakemake.output[0]


def get_smiles_from_chebi(chebi_id):
    numeric_id = chebi_id.replace("CHEBI:", "")
    url = f"https://www.ebi.ac.uk/chebi/backend/api/public/compound/{numeric_id}/"

    response = requests.get(url, timeout=30)
    response.raise_for_status()

    data = response.json()
    structure = data.get("default_structure")

    if structure and structure.get("smiles"):
        return structure["smiles"]

    return None


with open(fasta_path) as f:
    header = f.readline().strip()

match = re.match(r">[^|]+\|([^|]+)\|", header)
uniprot_id = match.group(1) if match else header.replace(">", "").split()[0]

url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
response = requests.get(url, timeout=60)
response.raise_for_status()

ligands = set()
current_ligand = None

for line in response.text.splitlines():
    m_name = re.search(r'/ligand="([^"]+)"', line)
    if m_name:
        current_ligand = m_name.group(1)

    m_id = re.search(r'/ligand_id="ChEBI:(CHEBI:\d+)"', line)
    if m_id and current_ligand:
        chebi_id = m_id.group(1)
        smiles = get_smiles_from_chebi(chebi_id)

        ligands.add((current_ligand, chebi_id, smiles))
        current_ligand = None

data = {
    "uniprot_id": uniprot_id,
    "ligands": [
        {
            "name": name,
            "chebi_id": chebi_id,
            "smiles": smiles,
        }
        for name, chebi_id, smiles in sorted(ligands)
    ],
}

with open(output_yaml, "w") as f:
    yaml.safe_dump(data, f, sort_keys=False)
