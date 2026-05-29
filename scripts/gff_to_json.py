#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import argparse
import math


def get_color_for_feature(feature_dict, colors, pathogen, subtype=None, segment=None):
    """
    Get the appropriate color for a feature based on the new colors.json structure.
    
    Colors structure:
    - Non-influenza: { pathogen: { gene_alias: [color_list] } }
    - Influenza: { pathogen: { subtype: { segment: [color_list] } } }
    
    For influenza segments with multiple colors, colors are assigned sequentially
    to CDS and mature_protein_region features within that segment.
    """
    if colors is None:
        return None

    pathogen_colors = colors.get(pathogen, None)
    if pathogen_colors is None:
        return None

    if pathogen in ["Influenza_A", "Influenza_B"]:
        if subtype is None:
            return None
        subtype_colors = pathogen_colors.get(subtype, None)
        if subtype_colors is None:
            return None
        if segment is None:
            return None
        segment_colors = subtype_colors.get(segment.lower(), None)
        if segment_colors is None:
            return None
        # Return the list of colors for this segment (caller will pick by index)
        return segment_colors
    else:
        # Non-influenza: lookup by alias/gene name
        if segment is not None:
            color_list = pathogen_colors.get(segment, None)
            if color_list and len(color_list) > 0:
                return color_list
        return None


def get_aliases_for_product(product, aliases, pathogen, subtype=None, segment=None):
    """
    Look up aliases for a product name based on the aliases.json structure.
    
    Aliases structure:
    - Non-influenza: { pathogen: { product_name: [aliases] } }
    - Influenza: { pathogen: { subtype: { segment: { product_name: [aliases] } } } }
    
    Returns the full list of aliases, or None if not found.
    """
    if aliases is None or product is None:
        return None

    pathogen_aliases = aliases.get(pathogen, None)
    if pathogen_aliases is None:
        return None

    if pathogen in ["Influenza A", "Influenza B"]:
        if subtype is None:
            return None
        subtype_aliases = pathogen_aliases.get(subtype, None)
        if subtype_aliases is None:
            return None
        # Search specified segment first
        if segment is not None:
            segment_aliases = subtype_aliases.get(segment, None)
            if segment_aliases is not None:
                alias_list = segment_aliases.get(product, None)
                if alias_list and len(alias_list) > 0:
                    return alias_list
        # Fallback: search all segments
        for seg_name, seg_aliases in subtype_aliases.items():
            if isinstance(seg_aliases, dict):
                alias_list = seg_aliases.get(product, None)
                if alias_list and len(alias_list) > 0:
                    return alias_list
        return None
    else:
        # Non-influenza: direct lookup
        alias_list = pathogen_aliases.get(product, None)
        if alias_list and len(alias_list) > 0:
            return alias_list
        return None


def add_protein_coordinates(feature_dict, pathogen="SARS-CoV-2"):
    """Add protein coordinates to the dictionary of features"""
    if feature_dict["type"] in ["mature_protein_region_of_CDS", "signal_peptide_region_of_CDS"]:
        if ":" in feature_dict.get("ID", ""):
            coord_part = feature_dict["ID"].split(":")[1]
            if ".." in coord_part:
                aa_start = int(coord_part.split("..")[0])
                aa_end = int(coord_part.split("..")[1])
                feature_dict["aa_start"] = aa_start
                feature_dict["aa_end"] = aa_end
    elif feature_dict["type"] == "CDS":
        aa_start = 1
        aa_end = math.floor((feature_dict["end"] - feature_dict["start"]) / 3)
        feature_dict["aa_start"] = aa_start
        feature_dict["aa_end"] = aa_end

    return feature_dict


def gff_to_json(gff_file_path, json_file_path, colors, aliases, pathogen, subtype=None, segment=None):
    """
    Convert GFF file to JSON format.
    
    Args:
        gff_file_path: Path to input GFF file
        json_file_path: Path to output JSON file
        colors: Parsed colors dictionary (new multi-pathogen format)
        aliases: Parsed aliases dictionary (new multi-pathogen format)
        pathogen: Pathogen name matching keys in colors/aliases JSON
        subtype: Subtype for influenza (e.g., H1N1, H3N2, Victoria)
        segment: Segment name for influenza (e.g., HA, NA, M2, NS)
    """
    # Map pathogen arg to the keys used in color and alias JSON files
    # Colors JSON uses underscores (Influenza_A), aliases uses spaces (Influenza A)
    color_pathogen_key = pathogen
    alias_pathogen_key = pathogen
    if pathogen == "Influenza_A":
        alias_pathogen_key = "Influenza A"
    elif pathogen == "Influenza_B":
        alias_pathogen_key = "Influenza B"
    elif pathogen == "SARS_CoV-2":
        alias_pathogen_key = "SARS-CoV-2"

    with open(gff_file_path, "r") as gff_file:
        features = {}
        id_list = []
        intergenic_dic = {}
        intergenic_dic["type"] = "INTERGENIC"
        intergenic_dic["color"] = "rgb(128,128,128)"

        species = None
        accession = None

        # Track color index for influenza segments with multiple colors
        segment_color_index = 0

        for line in gff_file:
            # Skip comment lines but extract metadata
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    accession = line.split(" ")[1].strip()
                elif line.startswith("##species"):
                    species = line.split(" ")[1].strip()
                continue

            fields = line.strip().split("\t")

            # Skip lines that don't have enough fields
            if len(fields) < 9:
                continue

            # Create feature dictionary
            feature_dict = {}
            feature_dict["type"] = fields[2]
            feature_dict["start"] = int(fields[3])
            feature_dict["end"] = int(fields[4])

            # Parse the attributes field
            attributes = {}
            for attribute in fields[8].split(";"):
                if attribute.strip() == "":
                    continue
                if "=" not in attribute:
                    continue
                key, value = attribute.split("=", 1)
                attributes[key] = value
            feature_dict.update(attributes)

            # Handle duplicate IDs
            if feature_dict.get("ID", "") in id_list:
                feature_dict["product"] = feature_dict.get("product", "unknown") + "-i"
            id_list.append(feature_dict.get("ID", ""))

            # Add protein coordinates
            feature_dict = add_protein_coordinates(feature_dict, pathogen=pathogen)

            # --- Resolve aliases (stored as a list) ---
            product = feature_dict.get("product", None)
            protein_aliases = get_aliases_for_product(
                product, aliases, alias_pathogen_key, subtype=subtype, segment=segment
            )
            if protein_aliases:
                feature_dict["protein_alias"] = protein_aliases

            # --- Resolve color ---
            if feature_dict["type"] == "CDS":
                # Determine the color key to look up
                color_lookup_key = None
                if protein_aliases:
                    # Use first alias as the color lookup key
                    color_lookup_key = protein_aliases[0]
                elif pathogen in ["Influenza_A", "Influenza_B"] and segment:
                    color_lookup_key = segment.lower()

                color_list = get_color_for_feature(
                    feature_dict, colors, color_pathogen_key,
                    subtype=subtype,
                    segment=segment if pathogen in ["Influenza_A", "Influenza_B"] else color_lookup_key
                )
                if color_list:
                    # Use first color for CDS
                    feature_dict["color"] = color_list[0]
                    segment_color_index = 1
                else:
                    feature_dict["color"] = "rgba(128, 128, 128, 1)"

            elif feature_dict["type"] in ["mature_protein_region_of_CDS", "signal_peptide_region_of_CDS"]:
                if pathogen in ["Influenza_A", "Influenza_B"]:
                    # For influenza, use sequential colors from the segment's color list
                    color_list = get_color_for_feature(
                        feature_dict, colors, color_pathogen_key,
                        subtype=subtype, segment=segment
                    )
                    if color_list and segment_color_index < len(color_list):
                        feature_dict["color"] = color_list[segment_color_index]
                        segment_color_index += 1
                    elif color_list:
                        feature_dict["color"] = color_list[-1]
                    else:
                        feature_dict["color"] = "rgba(128, 0, 128, 1)"
                else:
                    feature_dict["color"] = "rgba(128, 0, 128, 1)"

            elif feature_dict["type"] in ["stem_loop", "five_prime_UTR", "three_prime_UTR"]:
                feature_dict["color"] = "rgba(0, 0, 0, 1)"

            # Add species/accession to region features
            if feature_dict["type"] == "region":
                if species:
                    feature_dict["species"] = species
                if accession:
                    feature_dict["accession"] = accession

            # Print for debugging
            print(feature_dict)

            # --- Determine feature key ---
            feature_key = None

            if feature_dict["type"] == "CDS":
                if pathogen == "SARS_CoV-2":
                    if "protein_alias" in feature_dict:
                        feature_key = feature_dict["protein_alias"][0]
                    else:
                        if product and product == "ORF1a polyprotein":
                            continue  # Skip ORF1a polyprotein for SARS-CoV-2
                        elif product:
                            feature_key = product
                        else:
                            feature_key = feature_dict.get("ID", "unknown")
                else:
                    if "protein_alias" in feature_dict:
                        feature_key = feature_dict["protein_alias"][0]
                    elif product:
                        feature_key = product
                    elif "gene" in feature_dict:
                        feature_key = feature_dict["gene"]
                    else:
                        feature_key = feature_dict.get("ID", "unknown")

            elif feature_dict["type"] in ["mature_protein_region_of_CDS", "signal_peptide_region_of_CDS"]:
                if pathogen == "SARS_CoV-2":
                    if "protein_alias" in feature_dict and feature_dict["protein_alias"][0] not in features:
                        feature_key = feature_dict["protein_alias"][0]
                    elif "protein_alias" not in feature_dict and product and product in features:
                        feature_key = product
                    elif product:
                        feature_key = product
                    else:
                        feature_key = feature_dict.get("ID", "unknown")
                else:
                    if "protein_alias" in feature_dict:
                        feature_key = feature_dict["protein_alias"][0]
                    elif product:
                        feature_key = product
                    else:
                        feature_key = feature_dict.get("ID", "unknown")

            elif feature_dict["type"] == "gene":
                feature_key = feature_dict.get("ID", "unknown")

            elif "gbkey" in feature_dict:
                feature_key = feature_dict["gbkey"]
            else:
                feature_key = feature_dict.get("ID", "unknown")

            # Skip if feature_key is None
            if feature_key is None:
                continue

            # Add to features dictionary
            if pathogen == "SARS_CoV-2":
                # SARS-CoV-2: overwrite duplicates (original behavior)
                if feature_key not in features:
                    features[feature_key] = feature_dict
            else:
                # Other pathogens: keep all regions with unique keys
                if feature_key not in features:
                    features[feature_key] = feature_dict
                else:
                    suffix = 1
                    new_key = f"{feature_key}_{suffix}"
                    while new_key in features:
                        suffix += 1
                        new_key = f"{feature_key}_{suffix}"
                    features[new_key] = feature_dict

        # Add intergenic region
        features[intergenic_dic["type"]] = intergenic_dic

    # Write JSON output
    json_str = json.dumps(features, indent=4)
    with open(json_file_path, "w") as json_file:
        json_file.write(json_str)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Converts GFF to a dictionary of gene names and positions, "
        "and saves as JSON. Supports multiple pathogens including SARS-CoV-2, "
        "Measles, RSV-A, RSV-B, Influenza A, and Influenza B."
    )
    parser.add_argument("--gff_file", type=str, required=True,
                        help="Path to the GFF with genome annotation")
    parser.add_argument("--json_file", type=str, required=True,
                        help="JSON filename to save results to")
    parser.add_argument("--pathogen", type=str, required=True,
                        choices=["SARS_CoV-2", "Measles", "RSV-A", "RSV-B",
                                 "Influenza_A", "Influenza_B"],
                        help="Pathogen type (must match key in colors/aliases JSON)")
    parser.add_argument("--subtype", type=str, default=None,
                        help="Subtype for influenza (e.g., H1N1, H3N2, H5Nx, Victoria)")
    parser.add_argument("--segment", type=str, default=None,
                        help="Segment name for influenza (e.g., PB2, PB1, PA, HA, NP, NA, M2, NS)")
    parser.add_argument("--color_file", type=str, default=None,
                        help="JSON file containing color codes (new multi-pathogen format)")
    parser.add_argument("--alias_file", type=str, default=None,
                        help="JSON file containing alias mappings (new multi-pathogen format)")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()
    gff_file_path = args.gff_file
    json_file_path = args.json_file
    pathogen = args.pathogen
    subtype = args.subtype
    segment = args.segment

    # Validate influenza requires subtype
    if pathogen in ["Influenza_A", "Influenza_B"] and subtype is None:
        print(f"Error: --subtype is required for {pathogen}")
        exit(1)

    # Validate influenza requires segment
    if pathogen in ["Influenza_A", "Influenza_B"] and segment is None:
        print(f"Error: --segment is required for {pathogen}")
        exit(1)

    # Load colors
    colors = None
    if args.color_file:
        with open(args.color_file) as f:
            colors = json.load(f)

    # Load aliases
    aliases = None
    if args.alias_file:
        with open(args.alias_file) as f:
            aliases = json.load(f)

    # Check if the input file exists
    if not os.path.exists(gff_file_path):
        print(f"Error: Input file '{gff_file_path}' does not exist.")
        exit(1)

    # Convert GFF to JSON
    gff_to_json(gff_file_path, json_file_path, colors, aliases, pathogen,
                subtype=subtype, segment=segment)

    # Verify output
    if os.path.exists(json_file_path):
        print(f"Success: Output file '{json_file_path}' was created.")
    else:
        print(f"Error: Output file '{json_file_path}' was not created.")