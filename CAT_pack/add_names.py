#!/usr/bin/env python3

import argparse
import sys

import about
import check
import shared
import tax


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="CAT add_names",
        description="Add taxonomic names to CAT or BAT output files.",
        usage="CAT add_names -i FILE -o FILE -t DIR [options] [-h / --help]",
        add_help=False)
    
    required = parser.add_argument_group("Required arguments")
    shared.add_argument(
        required,
        "input_file",
        True,
        help_=("Path to input file. Can be classification or ORF2LCA output "
               "file from CAT or BAT."))
    shared.add_argument(required, "output_file", True)
    shared.add_argument(required, "taxonomy_folder", True)

    optional = parser.add_argument_group("Optional arguments")
    shared.add_argument(optional, "only_official", False)
    shared.add_argument(optional, "exclude_scores", False)
    shared.add_argument(optional, "force", False)
    shared.add_argument(optional, "quiet", False)
    shared.add_argument(optional, "help", False)

    (args, extra_args) = parser.parse_known_args()

    extra_args = [arg for (i, arg) in enumerate(extra_args) if
                  (i, arg) != (0, "add_names")]
    if len(extra_args) > 0:
        sys.exit("error: too much arguments supplied:\n{0}".format(
            "\n".join(extra_args)))

    # Add extra arguments.
    shared.expand_arguments(args)

    return args


def run():
    args = parse_arguments()

    message = "# CAT v{0}.".format(about.__version__)
    shared.give_user_feedback(
        message, args.log_file, args.quiet, show_time=False)

    errors = []

    errors.append(
        check.check_input_file(args.input_file, args.log_file, args.quiet))

    if not args.force:
        errors.append(
            check.check_output_file(
                args.output_file, args.log_file, args.quiet)
        )

    errors.append(
        check.check_in_and_output_file(
            args.input_file, args.output_file, args.log_file, args.quiet)
    )

    if True in errors:
        sys.exit(1)

    (taxid2parent,
            taxid2rank) = tax.import_nodes(
        args.nodes_dmp, args.log_file, args.quiet)
    taxid2name = tax.import_names(args.names_dmp, args.log_file, args.quiet)

    message = "Appending names..."
    shared.give_user_feedback(message, args.log_file, args.quiet)

    with open(args.input_file, "r") as f1:
        for line in f1:
            if line.startswith("#"):
                line = line.rstrip().split("\t")

                if "lineage" in line:
                    lineage_index = line.index("lineage")
                else:
                    message = ("{0} is not a supported classification file."
                               "".format(args.input_file))
                    shared.give_user_feedback(
                        message, args.log_file, args.quiet, error=True)

                    sys.exit(1)
                    
                try:
                    scores_index = line.index("lineage scores")
                except:
                    scores_index = None

                full_length = len(line)

                break
        else:
            message = ("{0} is not a supported classification file."
                       "".format(args.input_file))
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)
            
    with open(args.input_file, "r") as f1, open(args.output_file, "w") as outf1:
        for line in f1:
            line = line.rstrip()

            if line.startswith("#"):
                if args.only_official:
                    outf1.write("{0}\tsuperkingdom\tphylum\tclass\torder\t"
                                "family\tgenus\tspecies\n".format(line))
                else:
                    outf1.write("{0}\tfull lineage names\n".format(line))
                    
                continue
            
            line = line.split("\t")

            if len(line) != full_length:
                # Entry does not have a full annotation.
                outf1.write("{0}\n".format("\t".join(line)))

                continue

            if any([c.startswith("no taxid found") for c in line[2:4]]):
                # ORF has database hits but the accession number is not found
                # in the taxonomy files.
                outf1.write("{0}\n".format("\t".join(line)))

                continue
            
            lineage = line[lineage_index].split(";")

            if scores_index is not None and not args.exclude_scores:
                scores = line[scores_index].split(";")
            else:
                scores = None

            if args.only_official:
                names = tax.convert_to_official_names(
                    lineage, taxid2rank, taxid2name, scores
                )
            else:
                names = tax.convert_to_names(
                    lineage, taxid2rank, taxid2name, scores
                )

            outf1.write("{0}\t{1}\n".format("\t".join(line), "\t".join(names)))

    message = "Names written to {0}!".format(args.output_file)
    shared.give_user_feedback(message, args.log_file, args.quiet)

    return


if __name__ == "__main__":
    sys.exit("Run \'CAT add_names\' to add taxonomic names to CAT or BAT "
             "output files.")
