#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :

import sys, os
import networkx as nx
import argparse


def square_to_condensed(i, j, n):
    if j > i:
        tmp = j
        j = i
        i = tmp
    return int(n * i - ((i * (i + 1)) / 2) + j - 1 - i)


# See PopPUNK/trees.py
def build_rapidnj(rapidnj_exe, names, dists, out_prefix, threads=1):
    import io
    import subprocess

    try:
        from Bio import Phylo
    except ImportError as e:
        sys.stderr.write("Using rapidnj requires biopython")
        sys.exit(1)

    # generate phylip matrix
    phylip_name = f"{out_prefix}_core_distances.phylip"
    n_samples = len(names)

    with open(phylip_name, "w") as p_file:
        p_file.write(str(n_samples) + "\n")
        for i, ref_name in enumerate(names):
            p_file.write(ref_name)
            for j in range(n_samples):
                dist = dists[square_to_condensed(i, j, n_samples)] if i != j else 0.0
                p_file.write(f" {dist}")
            p_file.write("\n")

    # construct tree
    tree_filename = f"{out_prefix}.nwk"
    rapidnj_cmd = (
        f"{rapidnj_exe} {phylip_name} -n -i pd -o t -x {tree_filename} -c {threads}"
    )
    try:
        # run command
        subprocess.run(rapidnj_cmd, shell=True, check=True)

    # record errors
    except subprocess.CalledProcessError as e:
        sys.stderr.write(
            f'Could not run command "{rapidnj_cmd}"; returned code:{str(e.returncode)}\n'
        )
        sys.exit(1)

    # read tree
    tree = Phylo.read(tree_filename, "newick")

    # tidy unnecessary files
    os.remove(phylip_name)
    os.remove(tree_filename)

    # Midpoint root and format
    tree.root_at_midpoint()
    output_str = io.StringIO()
    Phylo.write(tree, output_str, "newick")
    output_str.seek(0)
    tree_string = output_str.getvalue()
    tree_string = tree_string.replace("'", "")

    return tree_string


## See PopPUNK/plot.py
def createMicroreact(prefix, pickle_loc, microreact_files, api_key=None):
    import pickle
    import requests
    import json
    from datetime import datetime

    microreact_api_new_url = "https://microreact.org/api/projects/create"
    description_string = "SKA distance run on " + datetime.now().strftime(
        "%Y-%b-%d %H:%M"
    )
    # Load example JSON to be modified
    with open(f"{pickle_loc}/microreact_example.pkl", "rb") as example_pickle:
        json_pickle = pickle.load(example_pickle)
    json_pickle["meta"]["name"] = description_string

    # Read data in
    with open(microreact_files[0]) as cluster_file:
        csv_string = cluster_file.read()
        json_pickle["files"]["data-file-1"]["blob"] = csv_string
    with open(microreact_files[1], "r") as dot_file:
        dot_string = dot_file.read()
        json_pickle["files"]["network-file-1"] = {
            "id": "network-file-1",
            "name": "network.dot",
            "format": "text/vnd.graphviz",
            "blob": dot_string,
        }
        json_pickle["networks"]["network-1"] = {
            "title": "Network",
            "file": "network-file-1",
            "nodeField": "id",
        }
    # Could add a tree later, for now this will be skipped
    if len(microreact_files) > 2:
        with open(microreact_files[2], "r") as tree_file:
            tree_string = tree_file.read()
            json_pickle["files"]["tree-file-1"]["blob"] = tree_string
    else:
        del json_pickle["files"]["tree-file-1"]
        del json_pickle["trees"]["tree-1"]

    with open(f"{prefix}.microreact", "w") as json_file:
        json.dump(json_pickle, json_file)

    url = None
    if api_key != None:
        headers = {
            "Content-type": "application/json; charset=UTF-8",
            "Access-Token": api_key,
        }
        r = requests.post(
            microreact_api_new_url, data=json.dumps(json_pickle), headers=headers
        )
        if not r.ok:
            if r.status_code == 400:
                sys.stderr.write(
                    "Microreact API call failed with response " + r.text + "\n"
                )
            else:
                sys.stderr.write(
                    "Microreact API call failed with unknown response code "
                    + str(r.status_code)
                    + "\n"
                )
        else:
            url = r.json()["url"]

    return url


def main():
    parser = argparse.ArgumentParser(
        prog="cluster_dists",
        description="Create clusters using results of `ska distance`",
        epilog="Requires networkx and pygraphviz",
    )
    parser.add_argument("distfile", help="Input distances (stdout from `ska distance`)")
    parser.add_argument("--output", help="Output prefix", default="ska_dist_clusters")
    parser.add_argument(
        "--snps", help="Maximum SNP distance to cluster", type=float, default=10
    )
    parser.add_argument(
        "--mismatches",
        help="Maximum k-mer mismatches to cluster",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--json-dir", help="Directory with 'microreact_example.pkl' file"
    )
    parser.add_argument(
        "--rapidnj", help="rapidnj executable (builds a tree from dists)"
    )
    parser.add_argument("--api-key", help="Microreact API key")

    args = parser.parse_args()

    build_tree = args.rapidnj != None
    if args.rapidnj != None:
        build_tree = True

    G = nx.Graph()
    samples = set()
    edges = list()
    dists = list()
    all_samples = False
    with open(args.distfile, "r") as distances:
        distances.readline()
        for line in distances:
            (sample1, sample2, snps, mismatches) = line.split("\t")
            # Don't have to keep adding samples once at the second row
            if not all_samples and sample1 in samples:
                all_samples = True
            else:
                samples.add(sample1)
                samples.add(sample2)
            if float(snps) <= args.snps and float(mismatches) <= args.mismatches:
                edges.append([sample1, sample2])
            if build_tree:
                dists.append(snps)

    # Create graph
    G.add_nodes_from(samples)
    G.add_edges_from(edges)
    microreact_files = list()

    # Write CSV of clusters
    with open(f"{args.output}.clusters.csv", "w") as out_csv:
        out_csv.write(",".join(["id", "Cluster__autocolour"]) + "\n")
        for idx, c in enumerate(
            sorted(nx.connected_components(G), key=len, reverse=True)
        ):
            for sample in c:
                out_csv.write(",".join([sample, str(idx + 1)]) + "\n")
    microreact_files.append(f"{args.output}.clusters.csv")

    # Layout and print (Microreact does the layout)
    # nx.nx_agraph.pygraphviz_layout(G, prog="neato")
    nx.nx_agraph.write_dot(G, f"{args.output}.graph.dot")
    microreact_files.append(f"{args.output}.graph.dot")

    if build_tree:
        tree = build_rapidnj(args.rapidnj, samples, dists, args.output)
        with open(f"{args.output}.njtree.nwk", "w") as tree_file:
            tree_file.write(tree)
        microreact_files.append(f"{args.output}.njtree.nwk")

    if args.json_dir == None:
        args.json_dir = "scripts/"
    url = createMicroreact(
        args.output,
        args.json_dir,
        microreact_files,
        args.api_key,
    )
    if url != None:
        sys.stderr.write("Microreact: " + url + "\n")
    else:
        sys.stderr.write("Provide --api-key to create microreact automatically\n")


if __name__ == "__main__":
    main()
