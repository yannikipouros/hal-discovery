This folder contains a fixed demo dataset used to validate hal-discovery.

It contains:
- 200 AF2 structures
- 10 radical halogenases (positive controls)
- 190 negative controls

To run:

hal-discovery mine \
  --pdb-dir Demo/Demo_dataset \
  --outdir Demo_output

Expected summary:
{
  "total_structures": 200,
  "total_2his_sites": 203,
  "proteins_with_2his_sites": 200,
  "sites_2his_aspglu": 192,
  "sites_putative_halogenase": 10
}
