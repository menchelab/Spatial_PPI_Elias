"""Command line interface"""
import argparse
from pathlib import Path
from typing import List

from io_helpers import *
from model import CellModel, Location
from config import load_config, build_constraint

def parse_list(s: str) -> List[str]:
    """convert CLI Argument into a list"""    
    if not s:
        return None
    items = [x.strip() for x in s.split(",")]
    items = [x for x in items if x]
    return items or None

def run() -> int:
    p = argparse.ArgumentParser(
        prog="cell_layout",
        description=(
            "Constrained 3D layout of a whole-cell PPI graph by locations/organelles. "
            "Inputs: PPI edge list + localization table + JSON config. "
            "Outputs: per-location 3D coordinates + optional Plotly graph."
        ),
    )

    p.add_argument("--ppi", required=True, help="Path to 2-column edge list TSV/CSV")
    p.add_argument("--localizations", required=True, help="Path to localization TSV")
    p.add_argument("--config", required=True, help="Path to JSON config (see examples)")
    p.add_argument("--outdir", required=True, help="Output directory")                  

    p.add_argument("--only", help="Comma-separated list of location names to include")      
    p.add_argument("--exclude", help="Comma-separated list of location names to exclude")   

    p.add_argument("--seed", type=int, help="Set seed for layout algorithm")                    
    p.add_argument("--iterations", type=int, help="Set iterations for layout")         

    p.add_argument("--plot", action="store_true", help="Interactive Plotly graph") 
    p.add_argument("--plot-title", default="3D cell", help="Plot title")          
    p.add_argument("--no-boundaries", action="store_false", help="Do not draw organelle boundaries") 
    p.add_argument("--no-edges", action="store_false", help="Do not draw edges")                     

    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load config and set defaults
    cfg = load_config(args.config)
    default_params = dict(cfg.default_params)
    
    if args.seed is not None:
        seed = args.seed
        print("set seed to ", seed)
    else:
        seed = default_params["seed"]
        print("used default seed ", seed)

    if args.iterations is not None:
        iterations = args.iterations
        print("set iterations to ", iterations)
    else:
        iterations = default_params["iterations"]
        print("used default iterations ", iterations)
    

    only = parse_list(args.only)
    exclude = set(parse_list(args.exclude) or [])
    
    # Load data
    G = load_ppi_graph(args.ppi)
    #print(G)
    loc_to_proteins, protein_to_locations, loc_df = load_localizations(args.localizations)
    init_locations_attribute(G)

    # Build model
    cell = CellModel(G)

    configured_names = []
    for spec in cfg.locations:
        if spec.name not in loc_to_proteins:
            print(f"Warning: {spec.name} is in the config but not found in the localization table.")
            continue

        if only is not None and spec.name not in only: #skip everything that is not in only
            continue
        if spec.name in exclude: #skip every excluded loaction
            continue

        raw = loc_to_proteins.get(spec.name, set())
        nodes = sorted(x for x in raw if isinstance(x, str) and x.strip() != "") # <- ignore NANs (occure because no loaction in locationfile)

        if not nodes:
            pass

        constraint = build_constraint(spec.constraint_type, spec.constraint_params)
        loc = Location(
            name=spec.name,
            node_ids=nodes,
            min_degree=spec.min_degree,
            center=spec.center,
            constraint=constraint,
            repulsion_strength=spec.repulsion_strength,
            seed=seed,
            iterations=iterations
        )
        cell.add_location(loc)
        configured_names.append(spec.name)

    if not configured_names:
        raise SystemExit("No locations selected. Check --only/--exclude and your config.")
    
    #print(cell.locations)

    # Run layout
    cell.assign_all_positions()

    #write positions in file
    out_path = write_positions(cell, outdir)
    print("Wrote:", out_path)

    # Summary
    cell.summary()
    
    #Plot
    if args.plot:
        cell.plot_all_locations_3d(title=args.plot_title, show_boundaries=args.no_boundaries, show_edges=args.no_edges)
        

    return 0

def main() -> None:
    raise SystemExit(run())


if __name__ == "__main__":
    main()


