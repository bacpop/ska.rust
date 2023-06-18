#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :

import sys, os
import networkx as nx
import argparse

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
    parser.add_argument("--api-key", help="Microreact API key")

    args = parser.parse_args()

    G = nx.Graph()
    samples = set()
    edges = list()
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

    # Create graph
    G.add_nodes_from(samples)
    G.add_edges_from(edges)

    # Layout and print (Microreact does the layout)
    # nx.nx_agraph.pygraphviz_layout(G, prog="neato")
    nx.nx_agraph.write_dot(G, f"{args.output}.graph.dot")

    # Write CSV of clusters
    with open(f"{args.output}.clusters.csv", "w") as out_csv:
        out_csv.write(",".join(["id", "Cluster__autocolour"]) + "\n")
        for idx, c in enumerate(
            sorted(nx.connected_components(G), key=len, reverse=True)
        ):
            for sample in c:
                out_csv.write(",".join([sample, str(idx + 1)]) + "\n")

    if args.json_dir == None:
        args.json_dir = "scripts/"
    url = createMicroreact(
        args.output,
        args.json_dir,
        [f"{args.output}.clusters.csv", f"{args.output}.graph.dot"],
        args.api_key,
    )
    if url != None:
        sys.stderr.write("Microreact: " + url + "\n")
    else:
        sys.stderr.write("Provide --api-key to create microreact automatically\n")


if __name__ == "__main__":
    main()
