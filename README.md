

# Spatial Representation of the Protein-Protein Interaction Network based on Associated Localization Data
### Build an organelle-aware, 3D visualization of a protein–protein interaction (PPI) network from a predefined data source and explore it in the MencheLab DataDiVR platform[1].


## Abstract
Using protein localization and PPI datasets, your team will (i) design a parameterized 3D “organelle space” for a single cell, (ii) map each protein to a compartment-specific (x, y, z) position using biologically meaningful placement rules, and (iii) render the full network in VR, where nodes are proteins positioned by their localization and edges represent physical interactions from the in-house PPI source. 

Your task is to define the relative organelle geometry and local placement logic (e.g., membrane-proximal vs. lumenal) or, where evidence is missing, assume organelle boundaries and areas within an absolute 3D space. Multi-localized proteins may be duplicated across compartments, subject to feasibility under the platform’s ~2 Mio. node constraints. If duplication would exceed limits, you will design an alternative encoding (e.g., toggles or weighted placement). Node annotations will serve the basis for interactive analysis and should be used to inform the network with additional data (e.g., organelle, evidence, confidence, function). Beyond visualization, the environment aims to support hypothesis generation e.g., assessing whether spatial proximity aligns with observed interaction patterns or pathway membership.

Deliverables include a VR “project” (DataDiVR-compatible file/folder structure), a fully reproducible Python pipeline (e.g. PEP 8 conform notebooks/scripts; Python ≥ 3.9.1), a documented git repository (template and example notebooks provided), and a final report (Introduction, Materials/Methods, Results, Discussion). 

MencheLab members will support data interpretation and domain questions throughout the semester. 
This project is suitable for 2-3 students. It will be supervised by the individual domain experts in the Menchelab 


## References 
(1) Pirch, S., Müller, F., Iofinova, E. et al. The VRNetzer platform enables interactive network analysis in Virtual Reality. Nat Commun 12, 2432 (2021). https://doi.org/10.1038/s41467-021-22570-w
