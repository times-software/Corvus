#!/usr/bin/env python3

import re
from collections import OrderedDict

##############################################################################
# CONFIGURATION
##############################################################################

KEYWORD_RE = re.compile(
    r"^[A-Za-z][A-Za-z0-9_.-]*$"
)

TYPE_MAP = {
    "S": "String",
    "I": "Integer",
    "F": "Float",
    "L": "Logical",
    "P": "Free-form"
}

SECTION_MAP = OrderedDict([

    ("Corvus Controls", [
        "target_list",
        "title",
        "write_input_only",
        "scratch",
        "usesaved",
        "multiprocessing_ncpu",
        "multiprocessing_level",
        "usehandlers"
    ]),

    ("General Methodology", [
        "method",
        "xc",
        "pspfiles",
        "scf_conv",
        "keep_symm",
        "mt_overlap",
        "molecule_vacuum_margin",
        "constant_volume",
        "nkpoints",
        "nqpoints",
        "pw_encut",
        "clusterradius"
    ]),

    ("Physical Properties", [
        "cif_input",
        "xyz_input",
        "xyz_snapshot",
        "vasp_",
        "cluster",
        "cell_",
        "species",
        "number_of_atoms",
        "absorbing_",
        "polarization",
        "charge",
        "atomic_charge",
        "multiplicity",
        "spectral_broadening",
        "fermi_shift",
        "spin_moment",
        "supercell_dimensions",
        "ismetal"
    ]),

    ("Abinit", [
        "abinit."
    ]),

    ("DMDW", [
        "dmdw."
    ]),

    ("NWChem", [
        "nwchem."
    ]),

    ("FEFF", [
        "feff.",
        "nuctemp",
        "debyetemp",
        "dmdw_nlanczos",
        "run_dmdw",
        "opcons."
    ]),

    ("ORCA", [
        "orca."
    ]),

    ("SIESTA", [
        "siesta.",
        "phsf."
    ]),

    ("cif2cell", [
        "cif2cell."
    ]),

    ("Loop Workflows", [
        "loop_"
    ]),

    ("Configurational Averaging", [
        "cfavg_"
    ]),

    ("Fitting", [
        "fit."
    ]),

    ("OCEAN", [
        "ocean."
    ])

])

##############################################################################
# HELPERS
##############################################################################

def categorize(keyword):

    for section, patterns in SECTION_MAP.items():

        for pattern in patterns:

            if keyword.startswith(pattern):
                return section

    return "Miscellaneous"


def latex_escape(text):

    if text is None:
        return ""

    text = str(text)

    replacements = {
        '\\': r'\textbackslash{}',
        '&': r'\&',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\textasciicircum{}',
    }

    for old, new in replacements.items():
        text = text.replace(old, new)

    return text


def make_label(keyword):

    return "kw:" + re.sub(
        r"[^A-Za-z0-9]+",
        "-",
        keyword
    )

##############################################################################
# SCHEMA ANALYSIS
##############################################################################

def parse_schema_info(schema_lines):

    rows = []

    for line in schema_lines:

        line = line.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        if "|X|" not in line:
            continue

        left, right = line.split("|X|", 1)

        types = []

        for token in right.split():

            if token in TYPE_MAP:
                types.append(token)

        rows.append({
            "default": left.strip(),
            "types": types
        })

    if not rows:

        return {
            "kind": "unknown",
            "default": "",
            "description": "Unknown"
        }

    repeated = any(
        "..." in line
        for line in schema_lines
    )

    #
    # Single-row schema
    #

    if len(rows) == 1:

        row = rows[0]

        if len(row["types"]) == 1:

            return {
                "kind": "scalar",
                "default": row["default"],
                "description":
                    TYPE_MAP[row["types"][0]]
            }

        return {
            "kind": "vector",
            "default": row["default"],
            "description":
                "Vector: "
                + ", ".join(
                    TYPE_MAP[t]
                    for t in row["types"]
                )
        }

    #
    # Repeating record
    #

    if repeated:

        signature = rows[0]["types"]

        return {
            "kind": "repeating-record",
            "default": "",
            "description":
                "Repeating record with columns: "
                + ", ".join(
                    TYPE_MAP[t]
                    for t in signature
                )
        }

    #
    # Matrix
    #

    signatures = {
        tuple(r["types"])
        for r in rows
    }

    if len(signatures) == 1:

        signature = rows[0]["types"]

        return {
            "kind": "matrix",
            "default": "",
            "description":
                f"{len(rows)} rows of "
                + ", ".join(
                    TYPE_MAP[t]
                    for t in signature
                )
        }

    #
    # Heterogeneous block
    #

    return {
        "kind": "block",
        "default": "",
        "description":
            "Multi-line heterogeneous block"
    }

##############################################################################
# PARSER
##############################################################################

def parse_parsnip(filename):

    with open(
        filename,
        "r",
        encoding="utf-8",
        errors="ignore"
    ) as f:

        lines = f.readlines()

    keywords = []

    i = 0

    while i < len(lines):

        line = lines[i].rstrip()

        #
        # Skip empty lines
        #

        if not line:
            i += 1
            continue

        #
        # Skip comment lines
        #

        if line.startswith("#"):
            i += 1
            continue

        #
        # Skip top-level documentation lines
        #

        if line.startswith("%"):
            i += 1
            continue

        #
        # Skip parser metadata
        #

        if ".inp_" in line:
            i += 1
            continue

        #
        # Skip schema content
        #

        if (
            line.startswith("{")
            or line.startswith("}")
            or line.startswith("|")
        ):
            i += 1
            continue

        keyword = line.strip()

        #
        # Validate keyword
        #

        if not KEYWORD_RE.fullmatch(keyword):
            i += 1
            continue

        ######################################################################
        # Documentation
        ######################################################################

        docs = []

        j = i + 1

        while j < len(lines):

            current = lines[j].rstrip()

            #
            # Documentation lines begin with %
            # Preserve as raw LaTeX.
            #

            if current.startswith("%"):

                docs.append(
                    current[1:].rstrip()
                )

                j += 1
                continue

            #
            # Beginning of schema
            #

            if current.strip().startswith("{"):
                break

            #
            # Any other non-empty line means
            # this keyword has ended.
            #

            if current.strip():
                break

            j += 1

        ######################################################################
        # Schema block
        ######################################################################

        schema_lines = []

        if (
            j < len(lines)
            and lines[j].strip().startswith("{")
        ):

            j += 1

            while j < len(lines):

                current = lines[j].rstrip()

                if "}" in current:
                    break

                schema_lines.append(current)

                j += 1

        ######################################################################
        # Analyze schema
        ######################################################################

        schema = parse_schema_info(
            schema_lines
        )

        ######################################################################
        # Deprecation detection
        ######################################################################

        deprecated = any(
            "deprecated" in doc.lower()
            for doc in docs
        )

        ######################################################################
        # Store keyword
        ######################################################################

        keywords.append({

            "keyword":
                keyword,

            "section":
                categorize(keyword),

            #
            # Raw LaTeX documentation lines
            #

            "documentation":
                docs,

            #
            # Original schema text
            #

            "schema":
                "\n".join(
                    schema_lines
                ),

            #
            # Schema analysis
            #

            "kind":
                schema["kind"],

            "default":
                schema["default"],

            "type_description":
                schema["description"],

            #
            # Metadata
            #

            "deprecated":
                deprecated,
        })

        #
        # Continue parsing after the schema
        #

        i = j

    return keywords


##############################################################################
# DOCUMENTATION WRITER
##############################################################################

def write_documentation(out, docs):

    if not docs:
        return

    for line in docs:

        #
        # Documentation lines are already valid LaTeX.
        # Do not escape anything.
        #

        out.write(line)
        out.write("\n")

    out.write("\n")

##############################################################################
# LATEX OUTPUT
##############################################################################

def write_latex(keywords, outfile):

    ##########################################################################
    # Group keywords by section
    ##########################################################################

    grouped = OrderedDict()

    for section in SECTION_MAP:
        grouped[section] = []

    grouped["Miscellaneous"] = []

    for kw in keywords:

        grouped.setdefault(
            kw["section"],
            []
        ).append(kw)

    ##########################################################################
    # Sort keywords within sections
    ##########################################################################

    for section in grouped:

        grouped[section] = sorted(
            grouped[section],
            key=lambda kw:
                kw["keyword"].lower()
        )

    ##########################################################################
    # Global alphabetical index
    ##########################################################################

    all_keywords = sorted(
        keywords,
        key=lambda kw:
            kw["keyword"].lower()
    )

    ##########################################################################
    # Open output file
    ##########################################################################

    with open(
        outfile,
        "w",
        encoding="utf-8"
    ) as out:

        ######################################################################
        # Preamble
        ######################################################################

        out.write(r"""
\documentclass[11pt]{article}

\usepackage[T1]{fontenc}
\usepackage{geometry}
\usepackage{hyperref}
\usepackage{longtable}
\usepackage{fancyvrb}
\usepackage{amsmath}
\usepackage{amssymb}

\geometry{margin=1in}

\title{Corvus Keyword Reference}
\author{Automatically Generated}
\date{\today}

\begin{document}

\maketitle
""")

        ######################################################################
        # Introduction
        ######################################################################

        out.write(r"""


\tableofcontents

\newpage

\section{Introduction}
\input{introduction}

\section{installation}
\input{installation}

\section{Running Corvus}
\input{running_corvus}

\section{The Corvus input file}

A Corvus input file is composed of a sequence of
keyword blocks.

Each keyword is followed by a brace-delimited
block containing the data associated with that
keyword.

Example:

\begin{Verbatim}
title
{
    Example Calculation
}
\end{Verbatim}

Supported schema types:

\begin{itemize}

\item S = String

\item I = Integer

\item F = Float

\item L = Logical

\item P = Free-form text

\end{itemize}

The documentation below lists every recognized
keyword, its expected schema, default value,
and documentation string.


""")

        ######################################################################
        # Section keyword index
        ######################################################################

        out.write(
            "\n\\section{Keyword Index by Section}\n"
        )

        for section, items in grouped.items():

            if not items:
                continue

            out.write(
                "\n\\subsection*{%s}\n"
                % latex_escape(section)
            )

            out.write(
                "\\begin{itemize}\n"
            )

            for kw in items:

                out.write(
                    "\\item "
                    "\\hyperref[%s]{%s}\n"
                    % (
                        make_label(
                            kw["keyword"]
                        ),

                        latex_escape(
                            kw["keyword"]
                        )
                    )
                )

            out.write(
                "\\end{itemize}\n"
            )

        ######################################################################
        # Alphabetical keyword index
        ######################################################################

        out.write(
            "\n\\section{Alphabetical Keyword Index}\n"
        )

        out.write(
            "\\begin{itemize}\n"
        )

        for kw in all_keywords:

            out.write(
                "\\item "
                "\\hyperref[%s]{%s}"
                " (%s)\n"
                % (
                    make_label(
                        kw["keyword"]
                    ),

                    latex_escape(
                        kw["keyword"]
                    ),

                    latex_escape(
                        kw["section"]
                    )
                )
            )

        out.write(
            "\\end{itemize}\n"
        )

        ######################################################################
        # Main reference
        ######################################################################

        for section, items in grouped.items():

            if not items:
                continue

            out.write(
                "\n\\section{%s}\n"
                % latex_escape(section)
            )

            for kw in items:

                ##################################################################
                # Keyword heading
                ##################################################################

                out.write(
                    "\n\\subsection{%s}\n"
                    "\\label{%s}\n"
                    % (
                        latex_escape(
                            kw["keyword"]
                        ),

                        make_label(
                            kw["keyword"]
                        )
                    )
                )

                ##################################################################
                # Metadata table
                ##################################################################

                out.write(
                    "\\begin{longtable}"
                    "{|p{1.5in}|p{4.5in}|}\n"
                    "\\hline\n"
                )

                rows = [

                    (
                        "Keyword",
                        kw["keyword"]
                    ),

                    (
                        "Kind",
                        kw["kind"]
                    ),

                    (
                        "Type",
                        kw["type_description"]
                    ),

                    (
                        "Default",

                        kw["default"]
                        if kw["default"]
                        else
                        "Required / No default"
                    ),

                    (
                        "Deprecated",

                        "Yes"
                        if kw["deprecated"]
                        else
                        "No"
                    ),
                ]

                for key, value in rows:

                    out.write(
                        "%s & %s \\\\\n"
                        "\\hline\n"
                        % (
                            latex_escape(key),
                            latex_escape(value)
                        )
                    )

                out.write(
                    "\\end{longtable}\n"
                )

                ##################################################################
                # Documentation
                ##################################################################

                write_documentation(
                    out,
                    kw["documentation"]
                )

                ##################################################################
                # Schema
                ##################################################################

                if kw["schema"]:

                    out.write(
                        "\\paragraph{Schema}\n"
                    )

                    out.write(
                        "\\begin{Verbatim}\n"
                    )

                    out.write(
                        kw["schema"]
                    )

                    out.write(
                        "\n\\end{Verbatim}\n"
                    )

        ######################################################################
        # Document end
        ######################################################################

        out.write(
            "\n\\end{document}\n"
        )

##############################################################################
# MAIN
##############################################################################

def main():

    keywords = parse_parsnip(
        "../corvutils/parsnip.corvus.config"
    )

    print(
        f"Parsed {len(keywords)} keywords"
    )

    write_latex(
        keywords,
        "Corvus_Manual.tex"
    )

    print(
        "Wrote Corvus_Manual.tex"
    )

if __name__ == "__main__":
    main()
