# Endothiapepsin to generate data for machine learning
######

# Aim
Zysiu requires molecular dynamics data.

# Chose of protein system

# Requirements for system

- ligand requirements:
    - small molecule
    - < 300 Da
    - rotatable bonds < 5
    - no "exotic" atoms like: S / F / Cl
- Protein
    - small
- Complex:
    - well definded binding pocket
    - clear interactions,
    - Hydrogen bonds
    - nice resolution
    - ligand should interact with ≥ 6 amino acids


# Data required

- see notion

# Protein preparation

1. Import 7IFH into maestro
2. All amino acids present
3. Remove GOL
4. Protein preparation:
   - pH = 4.6
   - default settings
   - remove water larger than 6 ang
  -> 2-7IFH-prep.maegz
5. Check quality
   1. Accept all first alternates in particular for Ile-125
   2. No water interactions with ligand
   -> 3-7IFH-clean.maegz

  # Complex description

  ## Key interactions


  ## Key amino acids

  - ASH-35
  - Asp-219         H-bond
  - Tyr-79:         pi-stacking
  - Thr-222
  - Leu-125         Hydrophobic
  - Gly-80          double H-bond
