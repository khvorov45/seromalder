from augur.distance import *


def vaccine_run(arg_alignment, distance_map_file, output, vaccine_strains):

    # Load sequences.
    alignments = load_alignments(arg_alignment, ["HA"])

    # Index sequences by node name and gene.
    sequences_by_node = defaultdict(dict)
    for gene, alignment in alignments.items():
        for record in alignment:
            sequences_by_node[record.name]["HA"] = str(record.seq)

    final_distances_by_node = {}
    distance_map_names = []

    attribute = "distance"

    # Load the given distance map.
    distance_map = read_distance_map(distance_map_file)
    distance_map_names.append(distance_map.get("name", distance_map_file))

    # Calculate distance between each sample and all other samples in a
    # previous season.
    distances_by_node = get_distances_to_vaccine_strains(
        sequences_by_node, distance_map, vaccine_strains
    )

    # Map distances to the requested attribute name.
    # Convert data like:
    # {
    #   "A/AbuDhabi/24/2017": 1
    # }
    # to data like:
    #
    # {
    #   "A/AbuDhabi/24/2017": {
    #     "ep": 1
    #   }
    # }
    #
    for node_name, values in distances_by_node.items():
        if node_name not in final_distances_by_node:
            final_distances_by_node[node_name] = {}

        final_distances_by_node[node_name][attribute] = values

    # Prepare params for export.
    params = {
        "attribute": attribute,
        "compare_to": "vaccine",
        "map_name": distance_map_names,
    }

    # Export distances to JSON.
    write_json({"params": params, "nodes": final_distances_by_node}, output)


def get_distances_to_vaccine_strains(sequences_by_node, distance_map, vaccine_strains):

    distances_by_node = {}

    # Calculate distance between each tip and all tips in previous seasons as
    # defined by the given latest date threshold.
    for vaccine_strain in vaccine_strains:  # tree.find_clades(terminal=True):

        vaccine_seq = sequences_by_node[vaccine_strain]
        if not vaccine_seq:
            print("vaccine name", vaccine_strain, "not found\n")
            continue

        # Distances between this node and other nodes are indexed by the other
        # node name.
        distances_by_node[vaccine_strain] = {}

        # Find all nodes that were sampled prior to the given latest date.
        for past_node_name in sequences_by_node.keys():

            past_node_seq = sequences_by_node[past_node_name]

            # Calculate distance between current node and the past node. Allow
            # comparison between any two nodes if an earliest or latest date is
            # not given.
            distances_by_node[vaccine_strain][
                past_node_name
            ] = get_distance_between_nodes(
                past_node_seq,
                vaccine_seq,
                distance_map,
            )

    return distances_by_node


vaccine_run(
    ["virus-distances/cdc-fasta-aligned.fasta"],
    "virus-distances/distance-map.json",
    "virus-distances/cdc-fasta-distances.json",
    [
        "EPI_ISL_165554|539576|A/Hong_Kong/4801/2014",
        "EPI_ISL_225834|780183|A/Singapore/INFIMH-16-0019/2016",
    ],
)
