# wasteflow-genome-configs

Reference genome configuration files for [WasteFlow 2.0](https://github.com/wasteflow/wasteflow) — a multi-pathogen wastewater genomics pipeline for probe-enriched sequencing data.

Downstream tools consume this repository by cloning it and resolving config paths via `index.json`.

---

## Repository structure

```
wasteflow-genome-configs/
├── index.json                          ← single source of truth; all paths resolve from here
├── schema/
│   ├── genome_config.schema.json       ← JSON Schema for non-segmented configs
│   └── segment_config.schema.json      ← JSON Schema for segmented configs
├── scripts/
│   ├── resolve_configs.py              ← library + CLI for path resolution
│   └── validate_index.py              ← CI validation (paths, orphans, duplicates)
└── viruses/
    ├── SARS-CoV-2/
    │   └── NC_045512.2/
    │       └── genome_config.json
    ├── RSV_A/
    │   └── PP109421.1/
    │       └── genome_config.json
    ├── RSV_B/
    │   └── OP975389.1/
    │       └── genome_config.json
    ├── Measles/
    │   └── NC_001498.1/
    │       └── genome_config.json
    ├── Influenza_A/
    │   ├── H1N1/
    │   │   ├── HA/CY121680.1/segment_config.json
    │   │   ├── NA/MW626056.1/segment_config.json
    │   │   └── ...
    │   ├── H3N2/
    │   │   └── ...
    │   └── H5Nx/
    │       └── ...
    └── Influenza_B/
        └── Victoria/
            └── ...
```

---

## Usage

### Clone the repository

```bash
git clone https://github.com/<org>/wasteflow-genome-configs.git
```

### Resolve config paths (Python library)

```python
from scripts.resolve_configs import GenomeConfigRegistry

registry = GenomeConfigRegistry("/path/to/wasteflow-genome-configs")

# Non-segmented virus
configs = registry.get_nonsegmented("SARS-CoV-2")
# → [("NC_045512.2", PosixPath("/path/to/.../genome_config.json"))]

# Single segment
configs = registry.get_segment("Influenza_A", "H1N1", "HA")
# → [("CY121680.1", PosixPath("/path/to/.../segment_config.json"))]

# All segments for a subtype
all_segs = registry.get_all_segments("Influenza_A", "H1N1")
# → {"HA": [...], "NA": [...], "PB1": [...], ...}

# List everything registered
registry.list_viruses()
```

### Resolve config paths (CLI)

```bash
# List all registered viruses
python scripts/resolve_configs.py --repo-root . --list

# Non-segmented
python scripts/resolve_configs.py --repo-root . --virus SARS-CoV-2

# All segments for a subtype
python scripts/resolve_configs.py --repo-root . --virus Influenza_A --subtype H1N1

# Specific segment
python scripts/resolve_configs.py --repo-root . --virus Influenza_A --subtype H1N1 --segment HA
```

---

## Adding a new virus or reference

### New non-segmented virus (e.g. Norovirus)

1. Create the config directory and file:
   ```
   viruses/Norovirus/MK032427.1/genome_config.json
   ```
2. Add an entry to `index.json` under `non_segmented`:
   ```json
   "Norovirus": {
     "display_name": "Norovirus GII",
     "references": [
       {
         "accession": "MK032427.1",
         "label": "GII.4 Sydney reference",
         "config": "viruses/Norovirus/MK032427.1/genome_config.json",
         "status": "active"
       }
     ]
   }
   ```
3. Run validation: `python scripts/validate_index.py --repo-root .`

### New reference for an existing virus (additional accession)

Append to the `references` array (non-segmented) or the segment accession array (segmented). Old accessions can be kept with `"status": "archived"` to retain history without exposing them to active pipelines.

### New Influenza subtype (e.g. H7N9)

1. Create segment directories under `viruses/Influenza_A/H7N9/<SEGMENT>/<ACCESSION>/`
2. Add a new subtype block under `segmented.Influenza_A.subtypes` in `index.json`

---

## Validation

The CI workflow (`.github/workflows/validate.yml`) runs on every push and PR to `main`. It checks:

- Every path in `index.json` exists on disk
- No config files on disk are absent from the index (orphan check)
- No duplicate accessions within the same scope

Run locally:

```bash
python scripts/validate_index.py --repo-root .
```

---

## Config file formats

Config files are generated from NCBI GFF annotations by WasteFlow's preprocessing scripts. See `schema/` for the JSON Schema definitions governing each format.

| File | Used for |
|---|---|
| `genome_config.json` | Non-segmented viruses (SARS-CoV-2, RSV, Measles) |
| `segment_config.json` | Individual segments of segmented viruses (Influenza A/B) |
