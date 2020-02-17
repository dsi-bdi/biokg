# -*- coding: utf-8 -*-


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


program_header = "= Building the [" + bcolors.OKBLUE + "biokg" + bcolors.ENDC + "] knowledge graph"
dwn_sym = "(" + bcolors.FAIL + "⤓" + bcolors.ENDC + ") "
done_sym = " (" + bcolors.OKGREEN + "✓" + bcolors.ENDC + ")"
fail_sym = " (" + bcolors.FAIL + "✘" + bcolors.ENDC + ")"
prc_sym = "(" + bcolors.OKBLUE + "⏣" + bcolors.ENDC + ") "
hsh_sym = " (" + bcolors.OKBLUE + "#" + bcolors.ENDC + ") "
inf_sym = "(" + bcolors.WARNING + bcolors.BOLD + "‣" + bcolors.ENDC + ") "


def print_line():
    """ print stdout line
    """
    print("------------------------------------------------")


def print_bold_line():
    """ print stdout line
    """
    print("================================================")


def print_section_header(header_txt):
    """

    Parameters
    ----------
    header_txt: str
        header text to print
    """
    print(">>>  %s ... " % header_txt)
    print_line()
