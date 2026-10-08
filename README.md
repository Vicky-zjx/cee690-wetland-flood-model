# Wetland Restoration and Flood Mitigation

**CEE690 — Systems Modeling for Climate Response, Duke University**

The project is physically organized into A, B, C and D. Each stage guide presents its final confirmed model, methods, results and file locations. All project descriptions are in English.

| Stage | Content | Guide |
|---|---|---|
| A | Initial eight-state model and illustrative synthetic baseline | [A — Initial model and baseline](A/README.md) |
| B | Fixed-gate wetland optimization and dynamic no-wetland comparison | [B — Current analysis](B/README.md) |
| C | Reserved; content not yet confirmed | [C](C/README.md) |
| D | Reserved; content not yet confirmed | [D](D/README.md) |

## Project folders

```text
Modeling/
├── README.md
├── .gitignore
├── A/
│   ├── README.md
│   ├── wetland_flood_model.py
│   ├── algebraic_plots.py
│   ├── model_notes.md
│   ├── requirements.txt
│   ├── plots/          A–J figures
│   ├── results/        Local baseline arrays and diagnostics
│   ├── reference/      Local course and model-reference material
│   └── work/           Shared official multivarious checkout
├── B/
│   ├── README.md
│   ├── *.py            Models, study drivers, reporting and checks
│   ├── *.json          Study configurations and source catalog
│   ├── data_sources_and_assumptions.md
│   ├── docs/           Detailed scientific reports
│   ├── plots/          Figures grouped by study
│   ├── tables/         Tables grouped by study
│   ├── results/        Local trajectories and numerical records
│   └── work/           Local observations, processing and provenance
├── C/
│   └── README.md
└── D/
    └── README.md
```

A and B are the physical stage folders, without assignment-number subfolders. Within B, names such as `wetland_only` and `wetland_counterfactual` identify different analyses. Wetland response labels A/A and B/B describe the two synthetic wetland classes; they do not subdivide the project stages.

## Execution and interpretation

Run all documented commands from this project root. Start with the [A environment setup](A/README.md#install-and-run). Both stages use Henri P. Gavin's official multivarious dependency, with the local checkout at `A/work/multivarious/`, revision `06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02`. No project license has been selected.

The current B study represents a hypothetical 0.66 km² installation driven by regional September 2022 event weather. At the finest tested timestep, its best-found wetland design lowers the reference river peak by **5.012 mm** under fixed gate operations. The separate no-wetland dynamic comparison is not another optimization. These uncalibrated simulations do not establish measured Milwaukee flood reductions or local engineering certification.

[.gitignore](.gitignore) keeps stage-local numerical results, working data, the external solver, course references, environments, caches and private material local. Saved configuration and provenance records retain their original path labels; [B/project_paths.py](B/project_paths.py) resolves those labels to the current physical locations. Saved numerical outputs are preserved. This organization is local; no GitHub publication is included.
