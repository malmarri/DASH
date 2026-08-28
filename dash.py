#!/usr/bin/env python3
"""
dash.py

(renamed from date_ancestral_segments.py 2026-08-28; same script)

Date a shared founder mutation from an *unphased* VCF: call the ancestral
segments, convert to genetic distance, and estimate its age -- this is the
whole point of the script, so --age-estimate runs by default whenever a
genetic map is given (pass --no-age-estimate to only do segment calling).

Along the way, this also builds the input file for the WEHI mutation-dating /
"Age of Mutation" app
(https://shiny.wehi.edu.au/rafehi.h/mutation-dating/) starting from an
*unphased* VCF.

The app's default input is a headerless, tab-separated, 3-column table:

    <individual ID>   <start bp of ancestral sharing>   <end bp of ancestral sharing>

one row per individual carrying the mutation.

How the segments are called here
--------------------------------
The documentation says ancestral segments are the region of "continuous
haplotype sharing between two or more individuals" as you move away from the
mutation, and that for individuals carrying two copies of the *same recessive*
mutation no phasing is needed: the segments are defined by continuous sharing
of identical homozygous markers.

This script implements exactly that, on unphased genotypes:

For every pair of carriers, walk outwards from the mutation. A marker is
concordant if BOTH individuals are homozygous there and homozygous for the
same allele; anything else (het in either, opposite homozygote) is discordant.
The pair stops sharing at the first run of --break-run consecutive discordant
markers. Each individual's segment is then the furthest point (left and right
independently) at which it still shares with at least --min-partners other
carrier(s) -- i.e. sharing between "two or more individuals".

This is not merely non-circular (an individual's own genotype never counts
toward its own boundary): checked directly against the documentation's own
worked example (Figure 1) -- the two longest-reaching individuals there stop
at the exact same marker, together, which is precisely max_j min(d_i, d_j),
this method's own formula, including its "longest capped at the
second-longest" consequence. An earlier consensus-haplotype method (each
individual walked against the majority homozygous allele across all carriers)
was tried and removed: it let an individual's own genotype vote toward the
"ancestral" haplotype it was then scored against, producing boundaries no
single other carrier actually reached -- a configuration the documentation's
own example never shows. See README.md "Method cross-check" for the retired
comparison and the reasoning.

Boundaries are reported at the *outermost concordant marker* (the last shaded
marker in Figure 1 of the documentation), which is what the app's cM
conversion expects. Use --boundary breakpoint to report the first discordant
marker instead (a more liberal, "up to the breakpoint" definition).

IMPORTANT (unphased data): this homozygosity-based approach is valid for
individuals homozygous for the same recessive mutation. Heterozygous carriers
(dominant / compound-het cases) cannot be handled this way -- the documentation
requires phased data for those. The script warns about any requested carrier
that is not homozygous for the ALT allele at the mutation position.

Extras
------
  * --report        per-individual QC (markers used, flanking markers, etc.)
  * --stats         median allele frequency + marker count for the chromosome,
                    i.e. two of the three "chance sharing" parameters of the
                    app's advanced options (the third, the genetic length of
                    the chromosome, must come from a genetic map).
  * --genetic-map   optionally also emit the advanced-option input: left- and
                    right-arm lengths in cM, one CSV row each, same individual
                    order, no row/column names.
  * --age-estimate  compute the age estimate itself (generations since the mutation,
                    both genealogy models, with/without the chance-sharing correction)
                    -- a stdlib port of the app's own Mutation_Age_Estimation.R, so the
                    whole pipeline (VCF -> segments -> cM arms -> age) runs in one step
                    without needing R or the Shiny app. On by default whenever
                    --genetic-map is given (that's the only other thing it needs);
                    pass --no-age-estimate to skip it.
  * --plot          publication-quality figure of the called haplotypes: one row
                    per carrier, one tick per marker coloured by the SAME
                    concordance test that called the segments, with the called
                    segment drawn underneath. Because it is generated from the
                    in-memory results it cannot drift from the numbers reported
                    above -- no re-reading of the VCF, no transcribed constants.

Requires only the Python 3 standard library, except --plot, which additionally
needs matplotlib.
"""

import argparse
import bisect
import csv
import gzip
import io
import math
import os
import statistics
import sys
from itertools import combinations

# genotype codes
HOM_REF, HET, HOM_ALT, MISSING = 0, 1, 2, -1


# --------------------------------------------------------------------------- #
# VCF reading
# --------------------------------------------------------------------------- #

def open_maybe_gzip(path):
    if path == "-":
        return sys.stdin
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def norm_chrom(c):
    c = str(c).strip()
    return c[3:] if c.lower().startswith("chr") else c


def load_mask(path, chrom):
    """Load a 3-column BED (chrom, start0, end) accessibility mask, keeping
    only intervals on `chrom` (matched with/without 'chr', like everywhere
    else in this script). Returns (starts, ends) as parallel sorted lists of
    0-based half-open interval bounds, merged if the input has overlaps.

    Equivalent in spirit to `bcftools view -R mask.bed`, done in-process so a
    single VCF pass can combine it with the other marker filters."""
    want = norm_chrom(chrom)
    ivs = []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or norm_chrom(f[0]) != want:
                continue
            ivs.append((int(f[1]), int(f[2])))
    ivs.sort()
    starts, ends = [], []
    for s, e in ivs:
        if starts and s <= ends[-1]:            # overlapping/adjacent -> merge
            ends[-1] = max(ends[-1], e)
        else:
            starts.append(s)
            ends.append(e)
    if not starts:
        sys.exit("ERROR: --mask has no intervals on chromosome %s" % chrom)
    return starts, ends


def in_mask(pos, starts, ends):
    """Is 1-based VCF position `pos` inside one of the (starts, ends) 0-based
    half-open intervals? O(log n) via bisect on the (non-overlapping,
    sorted) interval starts."""
    p0 = pos - 1
    i = bisect.bisect_right(starts, p0) - 1
    return i >= 0 and p0 < ends[i]


def het_ab_pvalue(ref_depth, alt_depth):
    """Exact two-sided binomial p-value for a heterozygous call's allele
    balance departing from 50:50 (allele-balance test only -- this looks at
    read counts, not strand; it will not catch strand-specific artefacts).

    Uses the standard "sum of outcomes no more probable than the one
    observed" definition, computed exactly with integer binomial
    coefficients (2**n cancels in the comparison so there's no float error
    until the final division). O(n) per call, n = ref_depth + alt_depth --
    only ever invoked on HET genotypes, so this is cheap in practice."""
    n = ref_depth + alt_depth
    if n == 0:
        return 1.0
    pk = math.comb(n, alt_depth)
    total = sum(c for i in range(n + 1) if (c := math.comb(n, i)) <= pk)
    # ldexp(total, -n) == total / 2**n but never builds 2.0**n as an
    # intermediate, which overflows for n >= 1024 (irrelevant at WGS depth,
    # but --min-dp is user-settable and amplicon/panel data can hit this).
    return min(1.0, math.ldexp(float(total), -n))


def gt_code(sample_field, gq_i=-1, dp_i=-1, ad_i=-1, min_gq=0.0, min_dp=0.0,
            het_ab_alpha=0.0, counter=None):
    """Return HOM_REF/HET/HOM_ALT/MISSING for a biallelic site, or None if the
    genotype references an allele index > 1 (marker should be dropped).

    A called genotype whose GQ/DP is present but below --min-gq/--min-dp is
    downgraded to MISSING, so a shaky call is treated as no information rather
    than as evidence against sharing. Genotypes with no GQ/DP at all are left
    alone -- that is how `bcftools merge -0` writes the hom-ref genotypes it
    fills in for sites absent from a single-sample VCF.

    A HET call whose AD ref:alt balance significantly departs from 50:50
    (--het-ab-alpha) is likewise downgraded to MISSING rather than treated as
    a hard discordance: this targets miscalled hets from reads that mapped
    confidently but to the wrong (paralogous/repetitive) locus, which GQ
    alone does not reliably flag. Deliberately depth-aware (an exact test,
    not a fixed ratio) so it does not over-flag real hets at low/moderate
    depth, where allele-balance sampling noise is large."""
    parts = sample_field.split(":")
    gt = parts[0]
    if not gt or gt[0] == ".":
        return MISSING
    alleles = gt.replace("|", "/").split("/")
    if any(a == "." or a == "" for a in alleles):
        return MISSING
    try:
        idx = [int(a) for a in alleles]
    except ValueError:
        return MISSING
    if any(i > 1 for i in idx):
        return None
    if len(idx) == 1:                      # haploid call (e.g. chrX males)
        code = HOM_REF if idx[0] == 0 else HOM_ALT
    elif idx[0] == idx[1]:
        code = HOM_REF if idx[0] == 0 else HOM_ALT
    else:
        code = HET
    for n, (i, thr) in enumerate(((gq_i, min_gq), (dp_i, min_dp))):
        if thr > 0 and 0 <= i < len(parts):
            v = parts[i].split(",")[0]
            if v not in (".", ""):
                try:
                    if float(v) < thr:
                        if counter is not None:
                            counter[n] += 1
                        return MISSING
                except ValueError:
                    pass
    if code == HET and het_ab_alpha > 0 and 0 <= ad_i < len(parts):
        adv = parts[ad_i].split(",")
        if len(adv) >= 2 and adv[0] not in (".", "") and adv[1] not in (".", ""):
            try:
                ref_d, alt_d = int(adv[0]), int(adv[1])
                if ref_d + alt_d > 0 and het_ab_pvalue(ref_d, alt_d) < het_ab_alpha:
                    if counter is not None:
                        counter[2] += 1
                    return MISSING
            except ValueError:
                pass
    if counter is not None:
        counter[3 if code == MISSING else 4] += 1
    return code


def read_vcf(path, chrom, args):
    """Stream the VCF, keep biallelic markers on `chrom`.

    Returns (samples, positions, codes, freqs) where
      positions[k] : bp of marker k (ascending)
      codes[k]     : list of genotype codes, one per sample, for marker k
      freqs[k]     : ALT allele frequency over all non-missing samples
    """
    want = norm_chrom(chrom)
    samples = None
    positions, codes, freqs = [], [], []
    n_seen = n_kept = 0
    d = args._filter_stats = {"multiallelic": 0, "filter": 0, "non_snv": 0,
                              "sym_or_N": 0, "allele_idx_gt1": 0, "no_call": 0,
                              "mask_excluded": 0,
                              "gq_set_missing": 0, "dp_set_missing": 0,
                              "het_ab_set_missing": 0,
                              "gt_already_missing": 0, "gt_kept": 0}

    mask = None
    if args.mask:
        mask_starts, mask_ends = load_mask(args.mask, chrom)
        mask = (mask_starts, mask_ends)
        sys.stderr.write("[info] --mask %s: %d accessible interval(s) on chromosome %s\n"
                         % (args.mask, len(mask_starts), chrom))

    fh = open_maybe_gzip(path)
    try:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if samples is None:
                sys.exit("ERROR: no #CHROM header line found in %s" % path)

            f = line.rstrip("\n").split("\t")
            if norm_chrom(f[0]) != want:
                continue
            n_seen += 1

            if mask is not None and not in_mask(int(f[1]), mask[0], mask[1]):
                d["mask_excluded"] += 1
                continue

            ref, alt, flt = f[3], f[4], f[6]
            if "," in alt:                                    # multiallelic
                d["multiallelic"] += 1
                continue
            if args.pass_only and flt not in ("PASS", ".", ""):
                d["filter"] += 1
                continue
            if args.snps_only and (len(ref) != 1 or len(alt) != 1):
                d["non_snv"] += 1
                continue
            if ref in ("N", "n") or alt in (".", "*"):
                d["sym_or_N"] += 1
                continue

            fmt = f[8].split(":")
            gq_i = fmt.index("GQ") if "GQ" in fmt else -1
            dp_i = fmt.index("DP") if "DP" in fmt else -1
            ad_i = fmt.index("AD") if "AD" in fmt else -1

            row = []
            bad = False
            n_alt = n_called = 0
            gt_counter = [0, 0, 0, 0, 0]
            for s in f[9:]:
                c = gt_code(s, gq_i, dp_i, ad_i, args.min_gq, args.min_dp,
                           args.het_ab_alpha, gt_counter)
                if c is None:
                    bad = True
                    break
                row.append(c)
                if c != MISSING:
                    n_called += 1
                    n_alt += 0 if c == HOM_REF else (1 if c == HET else 2)
            if bad:
                d["allele_idx_gt1"] += 1
                continue
            if n_called == 0:
                d["no_call"] += 1
                continue
            d["gq_set_missing"] += gt_counter[0]
            d["dp_set_missing"] += gt_counter[1]
            d["het_ab_set_missing"] += gt_counter[2]
            d["gt_already_missing"] += gt_counter[3]
            d["gt_kept"] += gt_counter[4]

            positions.append(int(f[1]))
            codes.append(row)
            freqs.append(n_alt / (2.0 * n_called))
            n_kept += 1
    finally:
        if fh is not sys.stdin:
            fh.close()

    if samples is None:
        sys.exit("ERROR: %s contains no VCF header" % path)
    if n_seen == 0:
        sys.exit("ERROR: no records found for chromosome %s (try --chrom with/without "
                 "the 'chr' prefix)" % chrom)
    if any(positions[i] > positions[i + 1] for i in range(len(positions) - 1)):
        order = sorted(range(len(positions)), key=lambda i: positions[i])
        positions = [positions[i] for i in order]
        codes = [codes[i] for i in order]
        freqs = [freqs[i] for i in order]

    sys.stderr.write("[info] chromosome %s: %d records read, %d markers kept after "
                     "filtering\n" % (chrom, n_seen, n_kept))
    sys.stderr.write("[info]   dropped: %d multiallelic, %d non-SNV, %d symbolic/N-ref, "
                     "%d FILTER, %d allele-index>1, %d no genotype call, %d outside --mask\n"
                     % (d["multiallelic"], d["non_snv"], d["sym_or_N"], d["filter"],
                        d["allele_idx_gt1"], d["no_call"], d["mask_excluded"]))
    if args.min_gq > 0 or args.min_dp > 0 or args.het_ab_alpha > 0:
        tot = (d["gt_kept"] + d["gq_set_missing"] + d["dp_set_missing"]
               + d["het_ab_set_missing"] + d["gt_already_missing"])
        sys.stderr.write("[info]   genotypes set to missing: %d by GQ<%g, %d by DP<%g, "
                         "%d by het allele-balance (p<%g) "
                         "(%.2f%% of %d genotypes; %d were already missing)\n"
                         % (d["gq_set_missing"], args.min_gq, d["dp_set_missing"],
                            args.min_dp, d["het_ab_set_missing"], args.het_ab_alpha,
                            100.0 * (d["gq_set_missing"] + d["dp_set_missing"]
                                     + d["het_ab_set_missing"]) / max(1, tot),
                            tot, d["gt_already_missing"]))
    return samples, positions, codes, freqs


# --------------------------------------------------------------------------- #
# marker selection
# --------------------------------------------------------------------------- #

def select_carriers(args, samples, positions, codes):
    """Resolve the carrier list and sanity-check genotypes at the mutation."""
    if args.samples_file:
        with open(args.samples_file) as fh:
            requested = [l.split()[0] for l in fh if l.strip() and not l.startswith("#")]
    elif args.samples:
        requested = [s for s in args.samples.replace(",", " ").split() if s]
    else:
        requested = None

    idx_at_mut = [k for k, p in enumerate(positions) if args.pos <= p <= args.mut_end]
    mut_codes = None
    if len(idx_at_mut) == 1:
        mut_codes = codes[idx_at_mut[0]]
    elif len(idx_at_mut) > 1:
        sys.stderr.write("[warn] %d markers overlap the mutation interval; carrier "
                         "auto-detection/checking uses the first one\n" % len(idx_at_mut))
        mut_codes = codes[idx_at_mut[0]]

    if requested is None:
        if mut_codes is None:
            sys.exit("ERROR: the mutation position is not a marker in the VCF, so carriers "
                     "cannot be auto-detected. Supply --samples or --samples-file.")
        requested = [s for s, c in zip(samples, mut_codes) if c == HOM_ALT]
        sys.stderr.write("[info] auto-detected %d carrier(s) homozygous ALT at %s:%d\n"
                         % (len(requested), args.chrom, args.pos))

    unknown = [s for s in requested if s not in samples]
    if unknown:
        sys.exit("ERROR: sample(s) not in VCF: %s" % ", ".join(unknown))

    col = {s: i for i, s in enumerate(samples)}
    carriers = [(s, col[s]) for s in requested]

    if mut_codes is not None:
        label = {HOM_REF: "hom REF", HET: "heterozygous", HOM_ALT: "hom ALT", MISSING: "missing"}
        odd = [(s, label[mut_codes[i]]) for s, i in carriers if mut_codes[i] != HOM_ALT]
        if odd:
            sys.stderr.write(
                "[warn] %d carrier(s) are not homozygous ALT at the mutation position: %s\n"
                "       Homozygosity-based segment calling assumes two copies of the same\n"
                "       recessive allele; het/compound-het carriers need phased data.\n"
                % (len(odd), ", ".join("%s=%s" % t for t in odd)))
    else:
        sys.stderr.write("[warn] the mutation position is not a marker in the VCF; "
                         "carrier genotypes were not checked\n")

    if len(carriers) < 2:
        sys.exit("ERROR: at least 2 carriers are required to call shared segments "
                 "(got %d)" % len(carriers))
    return carriers


def build_matrix(carriers, positions, codes, freqs, args):
    """Keep informative markers and return (pos, mat) with mat[k][i] the code of
    carrier i at marker k, plus the left/right walk orders."""
    n = len(carriers)
    cols = [i for _, i in carriers]
    max_missing = args.max_missing * n
    keep_pos, mat = [], []
    n_mut = n_maf = n_miss = 0
    for k, p in enumerate(positions):
        if args.pos <= p <= args.mut_end:      # the mutation itself is not evidence
            n_mut += 1
            continue
        if freqs[k] < args.min_maf or (1.0 - freqs[k]) < args.min_maf:
            n_maf += 1
            continue
        row = [codes[k][i] for i in cols]
        if sum(1 for c in row if c == MISSING) > max_missing:
            n_miss += 1
            continue
        keep_pos.append(p)
        mat.append(row)

    left = [k for k in range(len(keep_pos) - 1, -1, -1) if keep_pos[k] < args.pos]
    right = [k for k in range(len(keep_pos)) if keep_pos[k] > args.mut_end]
    sys.stderr.write("[info] %d informative markers around the mutation "
                     "(%d left, %d right)\n" % (len(keep_pos), len(left), len(right)))
    sys.stderr.write("[info]   dropped: %d in the mutation interval, %d below --min-maf, "
                     "%d with >%.0f%% of carriers missing\n"
                     % (n_mut, n_maf, n_miss, 100 * args.max_missing))
    if not left or not right:
        sys.stderr.write("[warn] no markers on one side of the mutation; that arm will "
                         "collapse to the mutation position\n")
    return keep_pos, mat, left, right


# --------------------------------------------------------------------------- #
# segment calling
# --------------------------------------------------------------------------- #

def walk(order, is_concordant, args, pos=None):
    """Walk markers outwards; return (last_concordant_k, first_break_k, n_used).

    Stops after --break-run consecutive discordant markers (or once
    --max-mismatch discordances have accumulated, if set).

    Discordant markers within --merge-mismatch-bp of the previous discordant
    marker count as ONE event: adjacent mismatches a few bases apart are one
    alignment/mapping artefact or one MNP, not independent evidence of
    recombination."""
    last_ok = None
    first_break = None
    last_bad_pos = None
    run = 0
    total = 0
    used = 0
    for k in order:
        state = is_concordant(k)
        if state is None:                       # uninformative (missing / non-hom pair)
            if args.missing_breaks:
                state = False
            else:
                continue
        used += 1
        if state:
            run = 0
            last_ok = k
            last_bad_pos = None      # a concordant marker ends any mismatch
                                      # cluster -- a later discordant marker
                                      # must not merge with one across it
        else:
            same_event = (args.merge_mismatch_bp > 0 and pos is not None
                          and last_bad_pos is not None
                          and abs(pos[k] - last_bad_pos) <= args.merge_mismatch_bp)
            if pos is not None:
                last_bad_pos = pos[k]
            if same_event:
                continue
            total += 1
            if first_break is None:
                first_break = k
            run += 1
            if run >= args.break_run or (args.max_mismatch >= 0 and total > args.max_mismatch):
                break
    return last_ok, first_break, used


def endpoint(last_ok, first_break, order, pos, fallback, args, side):
    """Translate a walk result into a base-pair coordinate."""
    if args.boundary == "breakpoint" and first_break is not None:
        k = first_break
    elif last_ok is not None:
        k = last_ok
    else:
        return fallback, None
    return pos[k], k


def call_pairwise(carriers, pos, mat, left, right, args):
    n = len(carriers)
    # per individual, per side: furthest bp reached with >= min_partners partners
    reach = {"L": [[] for _ in range(n)], "R": [[] for _ in range(n)]}
    detail = {"L": [[] for _ in range(n)], "R": [[] for _ in range(n)]}
    used_counts = [0] * n

    for i, j in combinations(range(n), 2):
        def conc(k, i=i, j=j):
            a, b = mat[k][i], mat[k][j]
            if a == MISSING or b == MISSING:
                return None
            if a == HET or b == HET:            # not homozygous -> sharing ends
                return False
            return a == b

        for side, order, fallback in (("L", left, args.pos), ("R", right, args.mut_end)):
            last_ok, first_break, used = walk(order, conc, args, pos)
            bp, k = endpoint(last_ok, first_break, order, pos, fallback, args, side)
            reach[side][i].append(bp)
            reach[side][j].append(bp)
            detail[side][i].append((bp, carriers[j][0]))
            detail[side][j].append((bp, carriers[i][0]))
            used_counts[i] = max(used_counts[i], used)
            used_counts[j] = max(used_counts[j], used)

    segments = []
    for i in range(n):
        m = max(1, args.min_partners)
        # sort (bp, partner_id) together so the reported partner is the one
        # whose reach actually set the m'th-ranked boundary, not always the
        # single most-distal partner (only differs from that when
        # --min-partners > 1; this project always runs with the default 1)
        l_detail_sorted = sorted(detail["L"][i], key=lambda t: t[0])
        r_detail_sorted = sorted(detail["R"][i], key=lambda t: -t[0])
        best_l = l_detail_sorted[m - 1] if len(l_detail_sorted) >= m else (args.pos, "-")
        best_r = r_detail_sorted[m - 1] if len(r_detail_sorted) >= m else (args.mut_end, "-")
        start, end = best_l[0], best_r[0]
        segments.append({
            "id": carriers[i][0], "start": start, "end": end,
            "left_partner": best_l[1], "right_partner": best_r[1],
            "markers_scanned": used_counts[i],
        })
    return segments


# --------------------------------------------------------------------------- #
# genetic map (optional advanced-option output)
# --------------------------------------------------------------------------- #

def load_map(path, chrom, fmt):
    """Return sorted [(bp, cM)] for `chrom`. Supports HapMap/1000G-style maps
    (chr pos rate cM, with or without header), PLINK .map (chr id cM bp) and
    a plain 'bp<sep>cM' two-column file."""
    want = norm_chrom(chrom)
    rows = []
    fh = open_maybe_gzip(path)
    try:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = line.replace(",", " ").split()
            low = [x.lower() for x in f]
            if any(x in ("position", "pos", "bp", "chr", "chromosome") for x in low) or \
               any("genetic_map" in x or x in ("cm", "cm_pos") for x in low):
                continue                                  # header line
            try:
                if fmt == "plink":
                    if norm_chrom(f[0]) != want:
                        continue
                    rows.append((int(float(f[3])), float(f[2])))
                elif fmt == "hapmap":
                    if len(f) >= 4:
                        if norm_chrom(f[0]) != want:
                            continue
                        rows.append((int(float(f[1])), float(f[3])))
                    else:                                  # pos rate cM
                        rows.append((int(float(f[0])), float(f[2])))
                else:                                      # two-column bp cM
                    rows.append((int(float(f[-2])), float(f[-1])))
            except (ValueError, IndexError):
                continue
    finally:
        if fh is not sys.stdin:
            fh.close()
    rows.sort()
    if len(rows) < 2:
        sys.exit("ERROR: genetic map %s yielded < 2 usable entries for chromosome %s "
                 "(try --map-format)" % (path, chrom))
    return rows


def interp_cm(rows, bp, warn):
    xs = [r[0] for r in rows]
    ys = [r[1] for r in rows]
    if bp <= xs[0] or bp >= xs[-1]:
        warn.add(bp)
        return ys[0] if bp <= xs[0] else ys[-1]
    j = bisect.bisect_left(xs, bp)
    if xs[j] == bp:
        return ys[j]
    x0, x1, y0, y1 = xs[j - 1], xs[j], ys[j - 1], ys[j]
    return y0 + (y1 - y0) * (bp - x0) / float(x1 - x0)


# --------------------------------------------------------------------------- #
# age estimation -- port of Mutation_Age_Estimation.R
# (Gandolfo, Bahlo & Speed 2014, Genetics 197:1315-1327;
#  github.com/bahlolab/DatingRareMutations), so the app's own age-estimation
# maths can be run in the same step as segment calling, without R.
# --------------------------------------------------------------------------- #

def reg_gammainc_lower(a, x):
    """Regularized lower incomplete gamma function P(a, x), i.e. the CDF of a
    Gamma(shape=a, scale=1) distribution at x. Standard series/continued-
    fraction algorithm (Numerical Recipes ch. 6.2); stdlib-only replacement
    for R's pgamma(). Accurate to ~1e-15."""
    if x < 0 or a <= 0:
        raise ValueError("bad a or x in reg_gammainc_lower")
    if x == 0:
        return 0.0
    gln = math.lgamma(a)
    if x < a + 1.0:                                # series representation
        ap = a
        summ = 1.0 / a
        delta = summ
        for _ in range(500):
            ap += 1.0
            delta *= x / ap
            summ += delta
            if abs(delta) < abs(summ) * 1e-15:
                break
        return summ * math.exp(-x + a * math.log(x) - gln)
    else:                                           # continued fraction, gives Q(a,x)
        tiny = 1e-300
        b = x + 1.0 - a
        c = 1.0 / tiny
        d = 1.0 / b
        h = d
        for i in range(1, 500):
            an = -i * (i - a)
            b += 2.0
            d = an * d + b
            if abs(d) < tiny:
                d = tiny
            c = b + an / c
            if abs(c) < tiny:
                c = tiny
            d = 1.0 / d
            delta = d * c
            h *= delta
            if abs(delta - 1.0) < 1e-15:
                break
        q = math.exp(-x + a * math.log(x) - gln) * h
        return 1.0 - q


def gamma_quantile(shape, prob, scale=1.0):
    """Inverse CDF of a Gamma(shape, scale) distribution at `prob` -- stdlib
    replacement for R's qgamma(p, shape, scale). Monotone bisection on the
    standardized (scale=1) variable, since reg_gammainc_lower is monotone
    increasing in x."""
    if prob <= 0.0:
        return 0.0
    if prob >= 1.0:
        return float("inf")
    lo, hi = 0.0, max(1.0, shape)
    while reg_gammainc_lower(shape, hi) < prob:
        hi *= 2.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if reg_gammainc_lower(shape, mid) < prob:
            lo = mid
        else:
            hi = mid
    return mid * scale


def estimate_age(l_lengths_cm, r_lengths_cm, confidence, chance_sharing_correction,
                  median_allele_frequency=None, markers_on_chromosome=None,
                  length_of_chromosome_cm=None):
    """Direct port of estimate_age() in Mutation_Age_Estimation.R (itself the
    body of that script) -- same variable names/order, so it can be checked
    line-by-line against the original. Returns generations-since-mutation
    point estimates and confidence intervals under both the 'independent' and
    'correlated' genealogy models (see the paper for the distinction)."""
    if chance_sharing_correction:
        e = 0.01
        p = median_allele_frequency ** 2 + (1 - median_allele_frequency) ** 2
        phi = (length_of_chromosome_cm / 100.0) / markers_on_chromosome
        loci = math.log(e) / math.log(p)
        cs_correction = loci * phi
    else:
        cs_correction = 0.0

    cc = confidence
    l = [x / 100.0 for x in l_lengths_cm]
    r = [x / 100.0 for x in r_lengths_cm]
    n = len(l)

    # --- 'independent' genealogy ---
    i_cs_correction = cs_correction if n >= 10 else 0.0
    length_correction = (sum(l) + sum(r) - 2 * (n - 1) * i_cs_correction) / (2 * n)
    sum_lengths = sum(l) + sum(r) + 2 * length_correction - 2 * (n - 1) * i_cs_correction
    b_c = (2 * n - 1) / (2.0 * n)
    i_tau_hat = (b_c * 2 * n) / sum_lengths
    g_l = gamma_quantile(2 * n, (1 - cc) / 2, 1.0 / (2 * n * b_c))
    g_u = gamma_quantile(2 * n, cc + (1 - cc) / 2, 1.0 / (2 * n * b_c))
    i_l, i_u = g_l * i_tau_hat, g_u * i_tau_hat

    # --- 'correlated' genealogy ---
    length_correction = (sum(l) + sum(r) - 2 * (n - 1) * cs_correction) / (2 * n)
    l2, r2 = list(l), list(r)
    l2[l2.index(max(l2))] += length_correction + cs_correction
    r2[r2.index(max(r2))] += length_correction + cs_correction
    lengths = [a + b - 2 * cs_correction for a, b in zip(l2, r2)]
    mean_len = statistics.mean(lengths)
    var_len = statistics.variance(lengths)          # sample variance (n-1), like R var()
    rho_hat = ((n * mean_len ** 2 - var_len * (1 + 2 * n))
               / (n * mean_len ** 2 + var_len * (n - 1)))
    n_star = n / (1 + (n - 1) * rho_hat)
    n_star = max(-n, min(n, n_star))
    b_c = (2 * n_star - 1) / (2.0 * n_star)
    c_tau_hat = (b_c * 2 * n) / sum(lengths)
    # NB: n.star may be reassigned below for strongly negative rho.hat, but
    # b.c (and hence the qgamma scale) is deliberately NOT recomputed after
    # that -- this matches the original R script exactly, quirk included.
    if rho_hat < -2.0 / (n - 1):
        n_star = n / (1 + (n - 1) * abs(rho_hat))
    if -2.0 / (n - 1) <= rho_hat < -1.0 / (n - 1):
        n_star = n
    g_l = gamma_quantile(2 * n_star, (1 - cc) / 2, 1.0 / (2 * n_star * b_c))
    g_u = gamma_quantile(2 * n_star, cc + (1 - cc) / 2, 1.0 / (2 * n_star * b_c))
    c_l, c_u = g_l * c_tau_hat, g_u * c_tau_hat

    return {"i_tau_hat": i_tau_hat, "i_l": i_l, "i_u": i_u,
            "c_tau_hat": c_tau_hat, "c_l": c_l, "c_u": c_u,
            "cs_correction": cs_correction, "n": n}


# --------------------------------------------------------------------------- #
# publication figure (--plot; the only part needing a non-stdlib library)
# --------------------------------------------------------------------------- #

FIG_SHARED = "#1b5e8c"      # dark blue  -- concordant
FIG_DIVERGENT = "#c0392b"   # red        -- opposite homozygote
FIG_HET = "#f0a30a"         # amber      -- heterozygous (lighter than the red so
                            #               the two stay separable under
                            #               red-green colour blindness)
FIG_MUT = "#6a1b9a"         # violet     -- mutation line; deliberately not red,
                            #               since discordant ticks are red and
                            #               both are thin vertical lines
FIG_BOUNDARY = "#1a1a1a"    # near-black -- called segment bar


def short_ids(ids, full=False):
    """Display labels: the text before the first '_' when that is still unique
    across carriers (typical of '<individual>_<run>_<batch>' sample names),
    otherwise the IDs unchanged."""
    if not full:
        cand = {s: s.split("_")[0] for s in ids}
        if len(set(cand.values())) == len(ids):
            return cand
    return {s: s for s in ids}


def make_figure(args, carriers, pos, mat, segments, map_rows=None,
                lefts=None, rights=None, age=None, color_test=None):
    """Draw the haplotype-sharing figure from the in-memory results.

    One row per carrier; one thin tick per informative marker, coloured by the
    SAME concordance test that called that individual's boundary: against the
    specific partner recorded in left_partner/right_partner, i.e. the partner
    whose comparison actually set that side's boundary (with --min-partners 1
    the segment end IS that partner's reach, so bars and colours agree by
    construction).

    color_test: optional callback (k, i, bp) -> (state, het) or None, for
    callers whose segment-calling method isn't this module's own pairwise walk
    (e.g. date_ancestral_segments_linkdatagen.py) -- lets the figure colour
    markers by that method's own concordance test instead of guessing at a
    lookalike. k is the index into `mat`/`pos` (NOT the windowed index used
    internally below); return None for an uninformative marker (skip it),
    or (state, het) where state is concordant/discordant and het marks a
    heterozygous discordance for its own colour. Ignored when left at its
    default of None, which draws the pairwise coloring described above.

    Discordant ticks are drawn ON TOP of concordant ones. Markers are typically
    denser than one per pixel column, so painting concordant last would hide a
    large share of the discordances -- and those are the rare, informative
    markers that terminate a segment.

    Isolated discordant ticks inside a called segment are expected, not error:
    --break-run and --merge-mismatch-bp tolerate a few without ending it.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from matplotlib.lines import Line2D
        import numpy as np
    except ImportError:
        sys.stderr.write("[warn] --plot needs matplotlib; figure skipped "
                         "(pip install --user matplotlib)\n")
        return

    ids = [c[0] for c in carriers]
    id_to_col = {s: i for i, s in enumerate(ids)}
    seg_by_id = {s["id"]: s for s in segments}
    label = short_ids(ids, args.plot_full_ids)

    # window: wide enough to show every called segment plus context
    if args.plot_window_bp:
        win = args.plot_window_bp
    else:
        reach = max(max(args.pos - s["start"], s["end"] - args.mut_end)
                    for s in segments)
        win = max(int(reach * 1.35), 1000)
    lo, hi = args.pos - win, args.mut_end + win
    keep = [k for k, bp in enumerate(pos) if lo <= bp <= hi]
    win_pos = [pos[k] for k in keep]
    win_mat = [mat[k] for k in keep]

    order = sorted(range(len(ids)),
                   key=lambda i: -(seg_by_id[ids[i]]["end"] - seg_by_id[ids[i]]["start"]))
    n = len(order)

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9.5,
                         "axes.linewidth": 0.8,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(9.5, 1.05 * n + 1.6))

    for row, i in enumerate(order):
        y = n - 1 - row
        seg = seg_by_id[ids[i]]
        lp = id_to_col.get(seg["left_partner"])
        rp = id_to_col.get(seg["right_partner"])
        sh, dv, ht = [], [], []
        for k, bp in enumerate(win_pos):
            if color_test is not None:                 # caller-supplied method
                r = color_test(keep[k], i, bp)
                if r is None:
                    continue
                state, het = r
            else:
                a = win_mat[k][i]
                j = lp if bp < args.pos else rp
                if j is None:
                    continue
                b = win_mat[k][j]
                if a == MISSING or b == MISSING:
                    continue
                het = (a == HET or b == HET)
                state = False if het else (a == b)
            xm = (bp - args.pos) / 1e6
            seg_xy = ((xm, y - 0.38), (xm, y + 0.38))
            (ht if het else (sh if state else dv)).append(seg_xy)

        # concordant underneath, discordances on top -- see docstring
        if sh:
            ax.add_collection(LineCollection(sh, colors=FIG_SHARED,
                                             linewidths=0.6, zorder=1.0))
        if dv:
            ax.add_collection(LineCollection(dv, colors=FIG_DIVERGENT,
                                             linewidths=0.7, zorder=1.5))
        if ht:
            ax.add_collection(LineCollection(ht, colors=FIG_HET,
                                             linewidths=0.9, zorder=2.0))

        # called segment: solid bar centred in the gap below this row's ticks
        xs = (seg["start"] - args.pos) / 1e6
        xe = (seg["end"] - args.pos) / 1e6
        bar_y = y - 0.5
        ax.plot([xs, xe], [bar_y, bar_y], color=FIG_BOUNDARY, linewidth=2.6,
                solid_capstyle="butt", zorder=3)
        for xb in (xs, xe):
            ax.plot([xb, xb], [bar_y - 0.06, bar_y + 0.06],
                    color=FIG_BOUNDARY, linewidth=1.2)

    ax.axvline(0, color=FIG_MUT, linewidth=1.1, linestyle="--", zorder=4)
    if args.mut_end > args.pos:                      # mutation spans >1 bp
        ax.axvspan(0, (args.mut_end - args.pos) / 1e6, color=FIG_MUT,
                   alpha=0.15, lw=0, zorder=0)
    ax.text(-0.026 * win / 1e6, n - 0.02, "mutation", color=FIG_MUT,
            ha="right", va="bottom", fontsize=8.5, fontstyle="italic")

    ax.set_xlim(-win / 1e6, (args.mut_end - args.pos + win) / 1e6)
    ax.set_ylim(-1.05, n + 0.35)
    ax.set_yticks(range(n))
    # Named partners (F11556, F15861, ...) dropped from the labels 2026-08-28: with
    # --min-partners 1 the named partner is just whichever carrier happened to reach
    # farthest marker-by-marker, which can hinge on a single missing genotype call
    # near a boundary (see README.md "Publication figure") -- confusing to name in
    # the figure itself. Tick coloring is still per-row, against that same partner;
    # the partner's identity is in AJH_chr3_report.final.tsv (left_partner/right_partner)
    # for anyone who needs it.
    ylabels = [label[seg_by_id[ids[order[n - 1 - r]]]["id"]] for r in range(n)]
    ax.set_yticklabels(ylabels, fontsize=8.5, linespacing=1.5)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Position relative to mutation (Mb)", labelpad=14)

    if map_rows:
        mbp = np.array([r[0] for r in map_rows], dtype=float)
        mcm = np.array([r[1] for r in map_rows], dtype=float)
        mut_cm = float(np.interp(args.pos, mbp, mcm))

        def bp_to_cm(x_mb):
            return np.interp(np.asarray(x_mb) * 1e6 + args.pos, mbp, mcm) - mut_cm

        def cm_to_bp(cm):
            return (np.interp(np.asarray(cm) + mut_cm, mcm, mbp) - args.pos) / 1e6

        secax = ax.secondary_xaxis("top", functions=(bp_to_cm, cm_to_bp))
        secax.set_xlabel("Genetic distance from mutation (cM)", labelpad=10)
        title_pad = 46
    else:
        title_pad = 14

    handles = [
        Line2D([0], [0], color=FIG_SHARED, lw=2.2,
               label="Concordant with %s (shared)" %
                     ("partner" if color_test is None else "majority")),
        Line2D([0], [0], color=FIG_DIVERGENT, lw=2.2, label="Discordant, opposite homozygote"),
        Line2D([0], [0], color=FIG_HET, lw=2.2, label="Discordant, heterozygous"),
        Line2D([0], [0], color=FIG_BOUNDARY, lw=2.6, label="Called shared segment"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=8.3, handlelength=1.6, columnspacing=1.3,
               bbox_to_anchor=(0.54, -0.07))

    lines = []
    if lefts is not None and rights is not None:
        lines.append("%.2f cM total shared" % (sum(lefts) + sum(rights)))
    if age is not None:
        # Correlated genealogy, point estimate + CI -- reverses the 2026-08-24 decision to
        # hide the CI (confirmed with the user 2026-08-28). That decision was because the
        # correlated branch sits 0.47% past the age model's n.star discontinuity and flips
        # which branch it lands on 76% of the time under a +/-1% arm-length perturbation
        # (check_model_stability.py / README "Correlated-model stability") -- still true,
        # so age["c_l"]/age["c_u"] below are not a stable/trustworthy interval, just what
        # the model reports on this exact run. Kept in the figure at the user's request;
        # the caveat lives in the README, not the plot.
        # "(correlated genealogy)" dropped 2026-08-28 to shorten the box -- it was
        # overlapping the mutation position's violet marker/dashed line. Still
        # correlated-genealogy's c_tau_hat/c_l/c_u under the hood; just not labelled
        # as such in the panel anymore.
        lines.append("Estimated Age: %.0f generations, 95%% CI (%.1f, %.1f)"
                     % (age["c_tau_hat"], age["c_l"], age["c_u"]))
    if lines:
        ax.text(0.985, 0.985, "\n".join(lines), transform=ax.transAxes,
                ha="right", va="top", fontsize=8, linespacing=1.6,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="#cccccc", alpha=0.92))

    title = args.plot_title or ("Shared ancestral haplotype at %s:%s"
                                % (args.chrom, format(args.pos, ",")))
    ax.set_title(title, fontsize=11.5, pad=title_pad)
    fig.tight_layout()

    if args.plot:
        outs = list(args.plot)
    else:                                    # --plot given without a value
        prefix = args.cm_prefix or os.path.splitext(args.out)[0]
        outs = [prefix + "_haplotype_sharing.pdf", prefix + "_haplotype_sharing.png"]
    for path in outs:
        fig.savefig(path, dpi=args.plot_dpi, bbox_inches="tight")
        sys.stderr.write("[info] wrote %s\n" % path)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def main():
    p = argparse.ArgumentParser(
        description="DASH: Dating Ancestral Shared Haplotypes from unphased VCFs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--vcf", required=True, help="input VCF (plain or bgzip/gzip, '-' for stdin)")
    p.add_argument("--chrom", required=True, help="chromosome of the mutation ('2' or 'chr2')")
    p.add_argument("--pos", required=True, type=int,
                   help="mutation position in bp (start position if it spans >1 bp)")
    p.add_argument("--mut-end", type=int, default=None,
                   help="end bp of the mutation, if it spans more than one base "
                        "[default: --pos]")
    p.add_argument("--samples", help="comma/space separated carrier IDs")
    p.add_argument("--samples-file", help="file with one carrier ID per line")
    p.add_argument("-o", "--out", default="ancestral_segments.txt",
                   help="output table: ID <tab> start <tab> end, no header")

    g = p.add_argument_group("marker filtering")
    g.add_argument("--mask", help="3-column BED (chrom, start0, end) of accessible regions "
                        "(e.g. the 1000G/HGDP strict accessibility mask); markers outside "
                        "it are dropped, equivalent to running `bcftools view -R mask.bed` "
                        "first but done in-process")
    g.add_argument("--pass-only", action="store_true", help="keep only FILTER=PASS records")
    g.add_argument("--snps-only", dest="snps_only", action="store_true", default=True,
                   help="keep only biallelic SNVs")
    g.add_argument("--allow-indels", dest="snps_only", action="store_false",
                   help="also use biallelic indels")
    g.add_argument("--min-maf", type=float, default=0.0,
                   help="drop markers with minor allele frequency below this (cohort MAF)")
    g.add_argument("--max-missing", type=float, default=0.2,
                   help="drop markers with more than this fraction of carriers missing")
    g.add_argument("--min-gq", type=float, default=0.0,
                   help="downgrade genotypes with GQ below this to missing (neutral); "
                        "genotypes carrying no GQ, e.g. hom-ref filled in by "
                        "'bcftools merge -0', are kept")
    g.add_argument("--min-dp", type=float, default=0.0,
                   help="downgrade genotypes with DP below this to missing (neutral)")
    g.add_argument("--het-ab-alpha", type=float, default=0.0,
                   help="downgrade a heterozygous call to missing (neutral) if its AD "
                        "ref:alt read balance significantly departs from 50:50 -- exact "
                        "two-sided binomial test, flagged when p < this; depth-aware, so "
                        "it does not over-flag real hets at low/moderate depth the way a "
                        "fixed ratio cutoff would; catches paralog/repeat mapping "
                        "artefacts that GQ alone can miss; 0 disables (default: 0.0)")

    g = p.add_argument_group("segment calling")
    g.add_argument("--break-run", type=int, default=2,
                   help="number of consecutive discordant markers that ends a segment "
                        "(2 tolerates isolated genotyping errors; 1 = no tolerance)")
    g.add_argument("--merge-mismatch-bp", type=int, default=0,
                   help="discordant markers within this many bp of the previous one "
                        "count as a single event (adjacent mismatches are usually one "
                        "artefact or MNP, not independent recombination evidence); "
                        "0 disables")
    g.add_argument("--max-mismatch", type=int, default=-1,
                   help="total discordant markers tolerated per arm (-1 = unlimited)")
    g.add_argument("--missing-breaks", action="store_true",
                   help="treat missing genotypes as discordant instead of skipping them")
    g.add_argument("--min-partners", type=int, default=1,
                   help="number of other carriers an individual must still share with "
                        "for the segment to continue")
    g.add_argument("--boundary", choices=("marker", "breakpoint"), default="marker",
                   help="report the outermost shared marker ('marker', as in the app's "
                        "Figure 1) or the first discordant marker ('breakpoint')")

    g = p.add_argument_group("extra output")
    g.add_argument("--report", nargs="?", const="-", default=None,
                   help="per-individual QC report (path, or omit value for stderr)")
    g.add_argument("--stats", action="store_true",
                   help="print chance-sharing parameters (median allele frequency and "
                        "number of markers on the chromosome)")
    g.add_argument("--genetic-map", help="genetic map for this chromosome; also writes "
                                         "left/right arm lengths in cM")
    g.add_argument("--map-format", choices=("hapmap", "plink", "two-column"), default="hapmap",
                   help="genetic map layout")
    g.add_argument("--cm-prefix", default=None,
                   help="prefix for the cM CSV files [default: --out without extension]")
    g.add_argument("--age-estimate", action=argparse.BooleanOptionalAction, default=True,
                   help="compute generations-since-mutation point estimates and "
                        "confidence intervals, under both the independent- and "
                        "correlated-genealogy models, with and without the chance-sharing "
                        "correction -- a direct port of the WEHI app's own "
                        "Mutation_Age_Estimation.R (Gandolfo, Bahlo & Speed 2014), run here "
                        "instead of via the Shiny app. On by default -- it's the point of "
                        "this script -- but needs --genetic-map to have anything to work "
                        "with, so it is silently skipped (not an error) if that isn't given; "
                        "pass --no-age-estimate to turn it off explicitly")
    g.add_argument("--confidence", type=float, default=0.95,
                   help="confidence coefficient for --age-estimate's intervals")
    g.add_argument("--chance-sharing-maf", type=float, default=None,
                   help="median population allele frequency for --age-estimate's "
                        "chance-sharing correction, e.g. from an external reference panel's "
                        "AF file. Overrides the value --stats computes from this VCF's own "
                        "samples, which is unreliable when there are few of them (e.g. a "
                        "VCF of carriers only -- allele frequency only takes values in "
                        "multiples of 1/(2*n_samples), and is upwardly biased by construction "
                        "wherever every sample is a carrier)")

    g = p.add_argument_group("figure")
    g.add_argument("--plot", nargs="*", default=None, metavar="FILE",
                   help="write a publication figure of the called haplotypes: one row "
                        "per carrier, one tick per marker coloured by the same "
                        "concordance test that called the segments, the called segment "
                        "drawn underneath, a cM axis on top if --genetic-map is given, "
                        "and the age estimate inset. Built from the in-memory results, "
                        "so it cannot disagree with the numbers printed above. Give a "
                        "one or more paths (each extension picks its own format, so "
                        "'--plot fig.pdf fig.png' renders once and writes both); with no "
                        "value, writes <prefix>_haplotype_sharing.pdf and .png. "
                        "Needs matplotlib")
    g.add_argument("--plot-window-bp", type=int, default=None,
                   help="half-width of the plotted window in bp [default: 1.35x the "
                        "longest called arm, so every segment is fully visible]")
    g.add_argument("--plot-title", default=None,
                   help="figure title [default: '<chrom>:<pos>']")
    g.add_argument("--plot-dpi", type=int, default=600,
                   help="raster resolution for --plot (ignored for PDF/SVG)")
    g.add_argument("--plot-full-ids", action="store_true",
                   help="label rows with full sample IDs instead of shortening them "
                        "to the text before the first underscore")

    args = p.parse_args()
    if args.mut_end is None:
        args.mut_end = args.pos
    if args.mut_end < args.pos:
        sys.exit("ERROR: --mut-end is before --pos")
    if args.age_estimate and not args.genetic_map:
        args.age_estimate = False       # nothing to compute it from; not an error since
                                        # --age-estimate defaults on and --genetic-map doesn't
        sys.stderr.write("[info] --age-estimate has nothing to work with without "
                         "--genetic-map -- skipped\n")

    samples, positions, codes, freqs = read_vcf(args.vcf, args.chrom, args)
    carriers = select_carriers(args, samples, positions, codes)
    pos, mat, left, right = build_matrix(carriers, positions, codes, freqs, args)

    segments = call_pairwise(carriers, pos, mat, left, right, args)

    with open(args.out, "w", newline="") as fh:
        for s in segments:
            fh.write("%s\t%d\t%d\n" % (s["id"], s["start"], s["end"]))
    sys.stderr.write("[info] wrote %d segments to %s\n" % (len(segments), args.out))

    if args.report is not None:
        rh = sys.stderr if args.report == "-" else open(args.report, "w")
        rh.write("id\tstart\tend\tlength_bp\tleft_arm_bp\tright_arm_bp\t"
                 "left_partner\tright_partner\tmarkers_scanned\n")
        for s in segments:
            rh.write("%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\n" % (
                s["id"], s["start"], s["end"], s["end"] - s["start"] + 1,
                args.pos - s["start"], s["end"] - args.mut_end,
                s["left_partner"], s["right_partner"], s["markers_scanned"]))
        if rh is not sys.stderr:
            rh.close()

    median_maf_all = None
    if args.stats or args.age_estimate:
        maf_all = sorted(min(f, 1.0 - f) for f in freqs)
        median_maf_all = statistics.median(maf_all)

    if args.stats:
        carrier_cols = set(i for _, i in carriers)
        non_carrier = [i for i in range(len(samples)) if i not in carrier_cols]
        line = ["", "Chance-sharing parameters for chromosome %s:" % args.chrom,
                "  Markers on chromosome     : %d (markers passing filters in this VCF)" % len(freqs),
                "  Median allele frequency   : %.4f (median MAF, all %d samples)"
                % (median_maf_all, len(samples))]
        if non_carrier:
            nc = []
            for k in range(len(positions)):
                n_alt = n_called = 0
                for i in non_carrier:
                    c = codes[k][i]
                    if c != MISSING:
                        n_called += 1
                        n_alt += 0 if c == HOM_REF else (1 if c == HET else 2)
                if n_called:
                    f = n_alt / (2.0 * n_called)
                    nc.append(min(f, 1.0 - f))
            if nc:
                line.append("  Median allele frequency   : %.4f (median MAF, %d non-carriers)"
                            % (statistics.median(sorted(nc)), len(non_carrier)))
        line += ["  Length of chromosome (cM) : not derivable from a VCF -- take it from the "
                 "genetic map",
                 "  NOTE: marker count/frequencies should reflect the marker set actually used "
                 "to call the segments.", ""]
        sys.stderr.write("\n".join(line) + "\n")

    map_rows = lefts = rights = age_nocs = None
    if args.genetic_map:
        map_rows = rows = load_map(args.genetic_map, args.chrom, args.map_format)
        warn = set()
        cm_mut_l = interp_cm(rows, args.pos, warn)
        cm_mut_r = interp_cm(rows, args.mut_end, warn)
        lefts, rights = [], []
        for s in segments:
            lefts.append(round(max(0.0, cm_mut_l - interp_cm(rows, s["start"], warn)), 6))
            rights.append(round(max(0.0, interp_cm(rows, s["end"], warn) - cm_mut_r), 6))
        prefix = args.cm_prefix or os.path.splitext(args.out)[0]
        for name, vals in (("left_arm_cM", lefts), ("right_arm_cM", rights)):
            path = "%s_%s.csv" % (prefix, name)
            with open(path, "w", newline="") as fh:
                csv.writer(fh).writerow(vals)
            sys.stderr.write("[info] wrote %s\n" % path)
        if warn:
            sys.stderr.write("[warn] %d position(s) fell outside the genetic map range and "
                             "were clamped to its ends\n" % len(warn))
        sys.stderr.write("\nAdvanced-option text input (same individual order as %s):\n"
                         "  Left arm lengths : %s\n  Right arm lengths: %s\n"
                         % (args.out, ", ".join("%g" % v for v in lefts),
                            ", ".join("%g" % v for v in rights)))
        chrom_cm = rows[-1][1] - rows[0][1]
        sys.stderr.write("  Genetic length of chromosome %s from this map: %.4f cM\n"
                         % (args.chrom, chrom_cm))

        if args.age_estimate:
            n = len(segments)
            sys.stderr.write("\nAge estimate (port of Mutation_Age_Estimation.R; %d "
                             "carriers, %.0f%% CI):\n" % (n, 100 * args.confidence))

            age_nocs = r_nocs = estimate_age(lefts, rights, args.confidence, False)
            sys.stderr.write(
                "  chance-sharing correction = FALSE (recommended, matches the app's default)\n"
                "    independent genealogy: %.1f generations, CI (%.1f, %.1f)\n"
                "    correlated  genealogy: %.1f generations, CI (%.1f, %.1f)\n"
                % (r_nocs["i_tau_hat"], r_nocs["i_l"], r_nocs["i_u"],
                   r_nocs["c_tau_hat"], r_nocs["c_l"], r_nocs["c_u"]))

            if args.chance_sharing_maf is not None:
                cs_maf = args.chance_sharing_maf
                maf_note = " (from --chance-sharing-maf)"
            else:
                cs_maf = median_maf_all
                maf_note = (" (from this VCF's own %d sample(s) -- unreliable with so few; "
                            "pass --chance-sharing-maf from an external reference panel if "
                            "available)" % len(samples))
            r_cs = estimate_age(lefts, rights, args.confidence, True,
                                median_allele_frequency=cs_maf,
                                markers_on_chromosome=len(freqs),
                                length_of_chromosome_cm=chrom_cm)
            gen_note = ("" if n >= 10 else
                       " (independent genealogy unaffected: the original R script's n<10 "
                       "guard only disables this correction for that branch, not the "
                       "correlated one)")
            sys.stderr.write(
                "  chance-sharing correction = TRUE (median MAF=%.4f%s, %d markers, "
                "%.4f cM chromosome)%s\n"
                "    independent genealogy: %.1f generations, CI (%.1f, %.1f)\n"
                "    correlated  genealogy: %.1f generations, CI (%.1f, %.1f)\n"
                % (cs_maf, maf_note, len(freqs), chrom_cm, gen_note,
                   r_cs["i_tau_hat"], r_cs["i_l"], r_cs["i_u"],
                   r_cs["c_tau_hat"], r_cs["c_l"], r_cs["c_u"]))

    if args.plot is not None:
        make_figure(args, carriers, pos, mat, segments, map_rows=map_rows,
                    lefts=lefts, rights=rights, age=age_nocs)


if __name__ == "__main__":
    main()
