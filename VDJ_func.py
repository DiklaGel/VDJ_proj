##############################################################################
#                         Functions for use with                             #
# TraCeR - a tool to reconstruct TCR sequences from single-cell RNA-seq data #
#                                                                            #
# Please see README and LICENCE for details of use and licence conditions.   #
# This software was written by Mike Stubbington (ms31@sanger.ac.uk) from the #
# Teichmann Lab, EMBL-EBI and WTSI (www.teichlab.org). Latest versions are   #
# available for download at www.github.com/teichlab/tracer.                  #
#                                                                            #
#      Copyright (c) 2015, 2016 EMBL - European Bioinformatics Institute     #
#      Copyright (c) 2016 Genome Research Ltd.                               #
#      Author: M.J.T. Stubbington ms31@sanger.ac.uk                          #
##############################################################################

from __future__ import print_function

import glob
import os
import re
import shutil
import subprocess
from collections import defaultdict, Counter
from time import sleep

import Levenshtein
import networkx as nx
import six
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from . import io

import copy


def process_chunk(chunk):
    store_VDJ_rearrangement_summary = False
    store_junction_details = False
    store_alignment_summary = False
    store_hit_table = False
    alignment_summary = []
    hit_table = []
    looking_for_end = False
    return_dict = defaultdict(list)
    for line_x in chunk:

        if store_VDJ_rearrangement_summary:
            VDJ_rearrangement_summary = line_x.split("\t")
            for i in VDJ_rearrangement_summary:
                return_dict['VDJ_rearrangement_summary'].append(i)
            store_VDJ_rearrangement_summary = False

        elif store_junction_details:
            junction_details = line_x.split("\t")
            for i in junction_details:
                return_dict["junction_details"].append(i)
            store_junction_details = False

        elif store_alignment_summary:
            if not line_x.startswith("#"):
                if line_x.startswith("Total"):
                    store_alignment_summary = False
                else:
                    return_dict['alignment_summary'].append(line_x)

        elif store_hit_table:
            if not looking_for_end:
                if not line_x.startswith("#"):
                    return_dict['hit_table'].append(line_x)
                    looking_for_end = True
            else:
                if line_x.startswith("#") or line_x.startswith("\n"):
                    store_hit_table = False
                else:
                    return_dict['hit_table'].append(line_x)

        elif line_x.startswith('# Query'):
            query_name = line_x.split(" ")[2]
            query_length = line_x.split(" ")[3]
            return_dict['query_length'] = int(query_length.split("=")[1])
            # return_dict['query_name'] = query_name

        elif line_x.startswith('# V-(D)-J rearrangement summary'):
            store_VDJ_rearrangement_summary = True

        elif line_x.startswith('# V-(D)-J junction details'):
            store_junction_details = True

        elif line_x.startswith('# Alignment summary'):
            store_alignment_summary = True

        elif line_x.startswith('# Hit table'):
            store_hit_table = True
    return (query_name, return_dict)



def is_rearrangement_productive(seq):
    # returns a tuple of three true/false values (productive, contains stop,
    # in-frame)
    seq_mod_3 = len(seq) % 3
    if seq_mod_3 == 0:
        in_frame = True
    else:
        in_frame = False

    seq = Seq(seq, IUPAC.unambiguous_dna)
    aa_seq = seq.translate()
    contains_stop = "*" in aa_seq

    if in_frame and not contains_stop:
        productive = True
    else:
        productive = False

    return (productive, contains_stop, in_frame)


def get_segment_name(name, pattern):
    match = pattern.search(name)
    number = match.group(1)
    if match.group(3):
        sub_number = match.group(3)
    else:
        sub_number = ""
    return (number)



def collapse_close_sequences(recombinants, locus):
    # pdb.set_trace()
    contig_names = [r.contig_name for r in recombinants]
    filtered_contig_names = [r.contig_name for r in recombinants]
    uncollapsible_contigs = []
    if len(recombinants) > 1:
        for i in range(len(recombinants) - 1):
            base_name = recombinants[i].contig_name
            base_seq = recombinants[i].imgt_reconstructed_seq
            base_V_segment = recombinants[i].best_VJ_names[0]
            base_J_segment = recombinants[i].best_VJ_names[1]

            base_id = recombinants[i].identifier
            base_junc = base_id.split("_")[1]
            base_e_value = float(recombinants[i].hit_table[0][-2])

            for j in range(i + 1, len(recombinants)):
                comp_name = recombinants[j].contig_name
                comp_seq = recombinants[j].imgt_reconstructed_seq
                comp_V_segment = recombinants[j].best_VJ_names[0]
                comp_J_segment = recombinants[j].best_VJ_names[1]
                comp_id = recombinants[j].identifier
                comp_junc = comp_id.split("_")[1]
                comp_e_value = float(recombinants[j].hit_table[0][-2])
                lev_dist = Levenshtein.distance(base_seq, comp_seq)
                # print("{}\t{}\t{}".format(base_id, comp_id, lev_dist))
                if lev_dist < 35 and not base_id == comp_id and base_name in filtered_contig_names \
                        and comp_name in filtered_contig_names:
                    # pdb.set_trace()
                    # define re pattern here to find TRAVx[DN] or TRDVx[DN]
                    # depending on locus
                    if locus == "TCR_A":
                        duplicate_pattern = re.compile(r"TRAV\d+[DN]")
                        segment_pattern = re.compile(
                            r"TRAV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    elif locus == "TCR_D":
                        duplicate_pattern = re.compile(r"DV\d+[DN]")
                        segment_pattern = re.compile(r"DV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    else:
                        uncollapsible_contigs.append(
                            "{}_vs_{}".format(base_name, comp_name))
                        attempt_collapse = False
                    if attempt_collapse and (
                        duplicate_pattern.search(
                            base_V_segment) or duplicate_pattern.search(
                            comp_V_segment)):
                        base_segment = get_segment_name(base_V_segment,
                                                        segment_pattern)
                        comp_segment = get_segment_name(comp_V_segment,
                                                        segment_pattern)
                        if base_segment == comp_segment:
                            # find alignment with lowest E value for V match
                            if base_e_value <= comp_e_value:
                                filtered_contig_names.remove(comp_name)
                            else:
                                filtered_contig_names.remove(base_name)
                        else:
                            uncollapsible_contigs.append(
                                "{}_vs_{}".format(base_name, comp_name))

                    else:
                        uncollapsible_contigs.append(
                            "{}_vs_{}".format(base_name, comp_name))

                elif lev_dist < 75 and not base_id == comp_id and base_name in filtered_contig_names \
                        and comp_name in filtered_contig_names:
                    if locus == "TCR_A":
                        duplicate_pattern = re.compile(r"TRAV\d+[DN]")
                        segment_pattern = re.compile(
                            r"TRAV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    elif locus == "TCR_D":
                        duplicate_pattern = re.compile(r"DV\d+[DN]")
                        segment_pattern = re.compile(r"DV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    else:
                        uncollapsible_contigs.append(
                            "{}_vs_{}".format(base_name, comp_name))
                        attempt_collapse = False
                    if attempt_collapse and (
                        duplicate_pattern.search(
                            base_V_segment) or duplicate_pattern.search(
                            comp_V_segment)):
                        base_segment = get_segment_name(base_V_segment,
                                                        segment_pattern)
                        comp_segment = get_segment_name(comp_V_segment,
                                                        segment_pattern)
                        if (base_segment == comp_segment) and (
                                base_junc == comp_junc) and (
                                base_J_segment == comp_J_segment):
                            # find alignment with lowest E value for V match
                            if base_e_value <= comp_e_value:
                                filtered_contig_names.remove(comp_name)
                            else:
                                filtered_contig_names.remove(base_name)
                        else:
                            uncollapsible_contigs.append(
                                "{}_vs_{}".format(base_name, comp_name))

                elif base_id == comp_id and base_name in filtered_contig_names and comp_name in filtered_contig_names:
                    if base_e_value <= comp_e_value:
                        filtered_contig_names.remove(comp_name)
                    else:
                        filtered_contig_names.remove(base_name)

    recombinants_to_delete = []

    for r in recombinants:
        if not r.contig_name in filtered_contig_names:
            recombinants_to_delete.append(r)

    [recombinants.remove(r) for r in recombinants_to_delete]

    return (recombinants)



def make_cell_network_from_dna(cells, keep_unlinked, shape, dot, neato,
                               receptor, loci,
                               network_colours):
    G = nx.MultiGraph()
    # initialise all cells as nodes

    if shape == 'circle':
        for cell in cells:
            G.add_node(cell, shape=shape,
                       label=cell.html_style_label_for_circles(receptor, loci,
                                                               network_colours),
                       sep=0.4, fontname="helvetica neue")
    else:
        for cell in cells:
            G.add_node(cell, shape=shape,
                       label=cell.html_style_label_dna(receptor, loci,
                                                       network_colours),
                       fontname="helvetica neue")
    # make edges:
    for i in range(len(cells)):
        current_cell = cells[i]
        comparison_cells = cells[i + 1:]

        for locus in loci:
            col = network_colours[receptor][locus][0]

            # current_identifiers = current_cell.getMainRecombinantIdentifiersForLocus(locus)
            for comparison_cell in comparison_cells:
                shared_identifiers = 0
                if current_cell.recombinants[receptor][locus] is not None:
                    for current_recombinant in \
                            current_cell.recombinants[receptor][locus]:
                        current_id_set = current_recombinant.all_poss_identifiers
                        if comparison_cell.recombinants[receptor][
                                locus] is not None:
                            for comparison_recombinant in \
                                    comparison_cell.recombinants[receptor][locus]:
                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                if len(current_id_set.intersection(
                                        comparison_id_set)) > 0:
                                    shared_identifiers += 1

                if shared_identifiers > 0:
                    width = shared_identifiers * 2
                    G.add_edge(current_cell, comparison_cell, locus,
                               penwidth=width, color=col,
                               weight=shared_identifiers)

    deg = G.degree()

    to_remove = [n for n in deg if deg[n] == 0]

    if len(to_remove) < len(G.nodes()):
        if not shape == 'circle':
            G.remove_nodes_from(to_remove)
            drawing_tool = [dot, '-Gsplines=true', '-Goverlap=false',
                            '-Gsep=0.4']
        else:
            drawing_tool = [dot, '-Gsplines=true', '-Goverlap=false']
    else:
        drawing_tool = [neato, '-Gsplines=true', '-Goverlap=false']

    bgcolors = ['#8dd3c720', '#ffffb320', '#bebada20', '#fb807220', '#80b1d320',
                '#fdb46220', '#b3de6920', '#fccde520',
                '#d9d9d920', '#bc80bd20', '#ccebc520', '#ffed6f20']
    component_counter = 0
    component_groups = list()
    j = 0
    components = nx.connected_components(G)

    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)
                G.node[cell]['style'] = 'filled'
                G.node[cell]['fillcolor'] = bgcolors[j]
                cell.bgcolor = bgcolors[j]

            if j < 11:
                j += 1
            else:
                component_counter += 1
                j = 0

        component_groups.append(members)

    return (G, drawing_tool, component_groups)


def draw_network_from_cells(cells, output_dir, output_format, dot, neato,
                            draw_graphs, receptor, loci, network_colours):
    cells = list(cells.values())
    network, draw_tool, component_groups = make_cell_network_from_dna(
        cells, False, "box", dot, neato, receptor, loci, network_colours)
    network_file = "{}/clonotype_network_with_identifiers.dot".format(
        output_dir)
    try:
        nx.write_dot(network, network_file)
    except AttributeError:
        import pydotplus
        nx.drawing.nx_pydot.write_dot(network, network_file)
    if draw_graphs:
        command = draw_tool + ['-o',
                               "{output_dir}/clonotype_network_with_identifiers.{output_format}".format(
                                   output_dir=output_dir,
                                   output_format=output_format), "-T",
                               output_format, network_file]
        subprocess.check_call(command)

    network, draw_tool, cgx = make_cell_network_from_dna(
        cells, False, "circle", dot, neato, receptor, loci, network_colours)
    network_file = "{}/clonotype_network_without_identifiers.dot".format(
        output_dir)
    try:
        nx.write_dot(network, network_file)
    except AttributeError:
        import pydotplus
        nx.drawing.nx_pydot.write_dot(network, network_file)
    if draw_graphs:
        command = draw_tool + ['-o',
                               "{output_dir}/clonotype_network_without_identifiers.{output_format}".format(
                                   output_dir=output_dir,
                                   output_format=output_format), "-T",
                               output_format, network_file]
        subprocess.check_call(command)
    return (component_groups)


def get_component_groups_sizes(cells, receptor, loci):
    cells = list(cells.values())
    G = nx.MultiGraph()
    # initialise all cells as nodes
    for cell in cells:
        G.add_node(cell)
    # make edges:
    for i in range(len(cells)):
        current_cell = cells[i]
        comparison_cells = cells[i + 1:]

        for locus in loci:

            for comparison_cell in comparison_cells:
                shared_identifiers = 0
                if current_cell.recombinants[receptor][locus] is not None:
                    for current_recombinant in \
                            current_cell.recombinants[receptor][locus]:
                        current_id_set = current_recombinant.all_poss_identifiers
                        if comparison_cell.recombinants[receptor][
                                locus] is not None:
                            for comparison_recombinant in \
                                    comparison_cell.recombinants[receptor][locus]:
                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                if len(current_id_set.intersection(
                                        comparison_id_set)) > 0:
                                    shared_identifiers += 1

                if shared_identifiers > 0:
                    width = shared_identifiers * 2
                    G.add_edge(current_cell, comparison_cell, locus,
                               penwidth=width, weight=shared_identifiers)

    deg = G.degree()

    to_remove = [n for n in deg if deg[n] == 0]

    # if len(to_remove) < len(G.nodes()):
    #    G.remove_nodes_from(to_remove)

    components = nx.connected_components(G)

    component_groups = list()

    singlets = []
    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)
            component_groups.append(members)
        else:
            for cell in component:
                singlets.append(cell.name)

    clonotype_size_counter = Counter([len(x) for x in component_groups])
    clonotype_size_counter.update({1: len(singlets)})

    clonotype_sizes = []
    max_size = max(list(clonotype_size_counter.keys()))
    if max_size < 5:
        for x in range(1, max_size + 1):
            clonotype_sizes.append(clonotype_size_counter[x])
        zero_padding = 5 - len(clonotype_sizes)
        clonotype_sizes = clonotype_sizes + [0] * zero_padding
    else:
        for x in range(1, max_size + 1):
            clonotype_sizes.append(clonotype_size_counter[x])

    return (clonotype_sizes)


def check_config_file(filename):
    if not os.path.isfile(filename):
        print()
        print("Couldn't find config file: {}".format(filename))
        print()
        exit(1)



def run_IgBlast(igblast, receptor, loci, output_dir, fasta, index_location,
                ig_seqtype, species,
                should_resume):
    print("##Running IgBLAST##")
    base_name = os.path.basename(fasta).split(".")[0]
    species_mapper = {
        'Mmus': 'mouse',
        'Hsap': 'human'
    }

    igblast_species = species_mapper[species]
    initial_locus_names = ["_".join([receptor, x]) for x in loci]
    locus_names = copy.copy(initial_locus_names)
    if should_resume:
        for locus in initial_locus_names:
            igblast_out = "{output_dir}/IgBLAST_output/{base_name}_{receptor}_{locus}.IgBLASTOut".format(
                output_dir=output_dir, base_name=base_name,
                receptor=receptor, locus=locus)
            if (os.path.isfile(igblast_out) and os.path.getsize(
                    igblast_out) > 0):
                locus_names.remove(locus)
                print(
                    "Resuming with existing IgBLAST output for {locus}".format(
                        locus=locus))

        if len(locus_names) == 0:
            return

    print("Performing IgBlast on {locus_names}".format(
        locus_names=locus_names))

    databases = {}
    for segment in ['V', 'D', 'J']:
        databases[segment] = "{}/{}_{}.fa".format(index_location, receptor,
                                                  segment)

    # Lines below suppress Igblast warning about not having an auxliary file.
    # Taken from
    # http://stackoverflow.com/questions/11269575/how-to-hide-output-of-subprocess-in-python-2-7
    DEVNULL = open(os.devnull, 'wb')

    for locus in locus_names:
        print("##{}##".format(locus))
        command = [igblast, '-germline_db_V', databases['V'],
                   '-germline_db_D', databases['D'],
                   '-germline_db_J', databases['J'], '-domain_system',
                   'imgt', '-organism', igblast_species,
                   '-ig_seqtype', ig_seqtype, '-show_translation',
                   '-num_alignments_V', '5',
                   '-num_alignments_D', '5', '-num_alignments_J', '5',
                   '-outfmt', '7', '-query', fasta]
        igblast_out = "{output_dir}/IgBLAST_output/{base_name}_{locus}.IgBLASTOut".format(
                output_dir=output_dir,
                base_name=base_name,
                locus=locus)
        with open(igblast_out, 'w') as out:
            # print(" ").join(pipes.quote(s) for s in command)
            subprocess.check_call(command, stdout=out, stderr=DEVNULL)

    DEVNULL.close()

