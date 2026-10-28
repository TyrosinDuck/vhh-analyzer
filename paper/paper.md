---
title: 'VHH tool box: a web application for single-domain antibody sequence and structure analysis'
tags:
  - Python
  - nanobody
  - VHH
  - single-domain antibody
  - antibody engineering
  - developability
  - phage display
  - bioinformatics
authors:
  - name: Jonathan Vaneyck
    orcid: 0000-0003-3095-7420
    affiliation: 1
affiliations:
  - name: Center for Protein Engineering (CIP), InBioS, University of Liège, Belgium
    index: 1
date: 28 October 2026
bibliography: paper.bib
---

# Summary

Single-domain antibodies — VHHs, or nanobodies — are the isolated variable
domains of the heavy-chain-only antibodies found in camelids
[@hamers-casterman1993]. At roughly 120 residues they are an order of magnitude
smaller than a conventional antibody, fold independently, and are produced in
bacteria, which has made them routine reagents in structural biology, imaging and
diagnostics [@muyldermans2013].

A laboratory that selects VHHs by phage display finishes a campaign holding
sequences, not answers. Which clones are genuinely different rather than the same
binder resampled? Which carry a deamidation site in a CDR, an unpaired cysteine,
or a hydrophobic patch that will cost yield? Which are close enough to human
germline to be worth humanising? Which of the genes that express poorly do so for
a reason visible in their codons? Each of these questions has an established
answer in the literature and a tool somewhere that computes it, and the
practical result is that the work is done by hand, in a spreadsheet, or in a
notebook that only its author can run.

`VHH tool box` answers them in one place, from a sequence. It classifies every
submitted protein and applies antibody-specific analyses only where they mean
something; annotates CDRs and frameworks in IMGT numbering [@lefranc2003] via
ANARCI [@dunbar2016]; scores the camelid hallmark positions and identity to the
closest human germline V gene; flags chemical liabilities and aggregation-prone
regions; measures solvent accessibility, disulfide bonds and antigen-contact
residues when a structural model is supplied; groups panning output into clonal
families by CDR3 identity at a threshold the user sets; back-translates to
*Escherichia coli*-optimised DNA and diagnoses existing genes for rare codons and
GC-rich stretches; extracts VHH sequences from Sanger chromatograms; and plans a
restriction cloning between two plasmids. It runs in a browser, requires no
installation, and stores nothing.

# Statement of need

The analyses above are individually well established and collectively
inaccessible. ANARCI numbers a domain but says nothing about developability;
developability predictors are built for therapeutic IgGs and assume a VL
interface that a VHH does not have; codon optimisers are generic and do not know
that the sequence is an antibody; and none of them read a chromatogram. A
researcher who wants all of it either installs five tools with incompatible
inputs or, far more commonly, does part of it and skips the rest.

`VHH tool box` is written for that researcher: someone who has just finished a
panning campaign, has between one and a hundred sequences, and needs to decide
what to order, express and keep. Its design follows from that use:

**Every number states how far to trust it.** The aggregation score is a stated
composite of charge, hydrophobicity, β-propensity and gatekeeper residues, not an
implementation of a parameterised predictor, and says so. The triage table
carries no composite score, deliberately: weighting a deamidation site against
three divergent framework positions depends on what the VHH is for, and a single
number would hide that choice behind a decimal point. Where a value is
theoretical, the interface says it is theoretical.

**Antibody assumptions are not applied to non-antibodies.** Each sequence is
classified, and a protein that is not a variable domain receives the generic
analyses only. Applied indiscriminately, the canonical-pair logic reports hen
lysozyme's eight ordinary cysteines as eight free thiols; such a result is worse
than no result, because it is confident.

**What a structure can measure is not inferred from sequence.** When a PDB,
mmCIF or AlphaFold model is supplied, SG–SG distances replace the canonical-pair
assumption, relative solvent accessibility distinguishes a buried liability motif
from an exposed one, hydrophobic patches are found by spatial contiguity rather
than by a sliding window over the sequence, and interface and epitope residues
are computed from the change in accessibility on binding — which sequence cannot
provide at all.

The software is in routine use at the Center for Protein Engineering, University
of Liège, where it replaced a set of private notebooks. It is distributed under
the MIT licence, documented in a user manual that states the limits of every
analysis, and covered by an automated test suite.

# Implementation

`VHH tool box` is written in Python and built on Biopython [@cock2009] for
sequence and structural operations, ANARCI [@dunbar2016] for IMGT numbering,
FAMSA [@deorowicz2016] for multiple alignment, `logomaker` [@tareen2020] for
sequence logos, and Streamlit for the interface. Relative solvent accessibility
is computed with the Shrake–Rupley algorithm and normalised by the residue maxima
of Tien et al. [@tien2013]; the codon adaptation index follows Sharp and Li
[@sharp1987] against a measured *E. coli* usage table. The synthetic-sequence
generator draws on a position-frequency matrix derived from PLAbDab-nano
[@plabdabnano].

The analysis logic is held in modules with no dependency on the interface, so
that the same functions can be called from a notebook, and the interface is a
thin layer of screens over them.

# Acknowledgements

The author thanks colleagues at the Center for Protein Engineering for testing
the application and reporting defects.

# References
