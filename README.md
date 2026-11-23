# PBI | Psychobiological Interaction

PBI is a research toolkit in development for PET-informed analysis of task fMRI. A reference PET map defines a spatial component, and a spatial regression estimates how strongly each subject's fMRI expresses that component over time.

The project focuses on two questions:

| Analysis | Outcome | Question |
| --- | --- | --- |
| Template-informed task response | A condition effect or contrast per subject and template | Does the expression of the PET-defined component change across task conditions? |
| Template-informed gPPI | A spatial map per subject, template, and contrast | Does coupling between that component and other brain locations differ across task conditions? |

Both analyses use fMRI-derived signals and a fixed PET spatial reference. The quantities being estimated are component response and functional coupling.

## Project layout

| Directory | Role |
| --- | --- |
| [matlab/](matlab/) | MATLAB/SPM implementation and reference workflow |
| [python/](python/) | Python implementation with no MATLAB/SPM runtime dependency |
| [resources/pet/](resources/pet/) | PET template catalog, provenance, and sources |
| [visualization/](visualization/) | Shared plots for quality control and results |
| [docs/](docs/methods.md) | Model specification and validation targets |

This is the initial development repository. The method specification and PET catalog are available now. The analysis implementations, surface support, and plotting modules are planned.

The first target is a volume workflow with a shared synthetic reference case for both languages. The surface extension will support a volume-defined seed with cortical targets. A seed estimated from the cortical surface is a separate analysis choice with different anatomical coverage.

## Related work

PBI builds on [REACT](https://doi.org/10.1016/j.neuroimage.2019.04.007) and [gPPI](https://doi.org/10.1016/j.neuroimage.2012.03.068). [Wong et al. (2022)](https://doi.org/10.1177/02698811221092509) and [Boucherie et al. (2023)](https://doi.org/10.1177/02698811231211154) describe related task applications.
