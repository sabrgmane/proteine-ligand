from pathlib import Path
import pandas as pd


input_csv = Path(snakemake.input[0])
output_csv = Path(snakemake.output[0])

df = pd.read_csv(input_csv)


score_columns = [
    "confidence_score",
    "ligand_iptm",
    "complex_iplddt",
    "L1_candidate_iptm",
    "L2_candidate_iptm",
    "L1_candidate_contacts_4A",
    "L2_candidate_contacts_4A",
    "L1_candidate_interface_residues",
    "L2_candidate_interface_residues",
    "L2_coordination_count_3A"
]

#normalisation
def normalize(series):
    plage = series.max() - series.min()
    if plage == 0:
        return pd.Series(1.0, ...)
    return (series - series.min()) / plage

def score_coordination_distance(series, optimal=2.1, tolerance=0.4):
    score = 1 - (series - optimal).abs() / tolerance
    return score.clip(lower=0, upper=1)


for col in score_columns:
    df[f"{col}_norm"] = normalize(df[col])

df["L2_mean_coordination_distance_3A_norm"] = score_coordination_distance(
    df["L2_mean_coordination_distance_3A"]
)

norm_cols = [f"{col}_norm" for col in score_columns] + [
    "L2_mean_coordination_distance_3A_norm"
]

#La moyenne n'est pas un score optimal, utilisée uniquement dans le cadre du test
df["ranking_score"] = df[norm_cols].mean(axis=1)

df = df.sort_values("ranking_score", ascending=False)
df.insert(0, "rank", range(1, len(df) + 1))

output_csv.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output_csv, index=False)
