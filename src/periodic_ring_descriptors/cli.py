"""Command-line interface for periodic ring descriptors."""

import argparse
import json
from pathlib import Path

import pandas as pd

from .descriptor_2d import compute_descriptor_2d
from .descriptor_3d import compute_descriptor_3d


def _write_result(result, json_out, rings_out=None):
    Path(json_out).write_text(json.dumps(result, indent=2), encoding="utf-8")
    if rings_out:
        pd.DataFrame(result["rings"]).to_csv(rings_out, index=False)


def main():
    parser = argparse.ArgumentParser(prog="ringdesc")
    sub = parser.add_subparsers(dest="dimension", required=True)

    p2 = sub.add_parser("2d", help="Periodic 2D ring descriptor")
    p2.add_argument("structure")
    p2.add_argument("--cutoff", type=float, default=1.895)
    p2.add_argument("--max-ring", type=int, default=12)
    p2.add_argument("--a", type=float, default=0.25)
    p2.add_argument("--json-out", default="descriptor_2d.json")
    p2.add_argument("--rings-out", default=None)

    p3 = sub.add_parser("3d", help="Periodic 3D edge-shortest ring descriptor")
    p3.add_argument("structure")
    p3.add_argument("--cutoff", type=float, default=1.895)
    p3.add_argument("--max-ring", type=int, default=20)
    p3.add_argument("--a", type=float, default=0.25)
    p3.add_argument("--image-range", type=int, default=1)
    p3.add_argument("--json-out", default="descriptor_3d.json")
    p3.add_argument("--rings-out", default=None)

    args = parser.parse_args()
    if args.dimension == "2d":
        result = compute_descriptor_2d(
            args.structure,
            cutoff=args.cutoff,
            max_ring_size=args.max_ring,
            a=args.a,
        )
    else:
        result = compute_descriptor_3d(
            args.structure,
            cutoff=args.cutoff,
            max_ring_size=args.max_ring,
            joukowsky_a=args.a,
            image_range=args.image_range,
        )
    _write_result(result, args.json_out, args.rings_out)
    print(f"Wrote {args.json_out}")


if __name__ == "__main__":
    main()