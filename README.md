# HipHap (previously Diplinator)

<img width="300" height="200" alt="hiphap logo" src="https://github.com/user-attachments/assets/794d2471-cdc8-46cb-b63a-f1c11fb5d51e" /> 

Diploid genome assemblies are now routinely available, but most read aligners were designed for haploid references. When reads are aligned to a diploid assembly, the aligner sees two nearly identical alignments to either haplotype, and thus reduces the mapping quality (MapQ) score to reflect this ambiguity. This can cause downstream tools to discard reads from easily mappable regions.
HipHap resolves this issue by aligning reads to each haplotype assembly separately and assigning each read to its best-supported haplotype. We also introduce a haplotype assignment quality score (HapQ) in HipHap to quantify confidence in the haplotype of origin of a read.

HipHap is implemented in Rust, and supports SAM, BAM, CRAM, and PAF formats


## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  - [Threads](#threads)
- [Example Workflow](#example-workflow)
  - [Diploid assembly alignment](#diploid-assembly-alignment)
  - [One file per haplotype (`-p`)](#one-file-per-haplotype--p)
  - [Comparing different reference genomes](#comparing-different-reference-genomes)
  - [CRAM input files](#cram-input-files)
- [Example PAF Usage](#example-paf-usage)
- [Weighted Alignment Scoring Mechanism](#weighted-alignment-scoring-mechanism)
- [Haplotype tag (HP)](#haplotype-tag-hp)
- [Losing-haplotype alignments (`--keep-loser`)](#losing-haplotype-alignments---keep-loser)
  - [Score threshold (`--loser-frac`)](#score-threshold---loser-frac)
  - [PAF](#paf)
- [Haplotype Assignment Quality (HapQ)](#haplotype-assignment-quality-hapq)
- [Citation](#citation)


## Installation

Recommended: Download precompiled binary:
```bash
# Todo
```

Build from source:
```bash
git clone https://github.com/jheinz27/hiphap.git
cd hiphap
cargo build --release
./target/release/hiphap
```

## Usage

```
HipHap: Choose the best alignment of a read to each haploid of a diploid assembly

Usage: hiphap [OPTIONS] <ASM1> <ASM2>

Arguments:
  <ASM1>  asm1 alignment file (sam/bam/cram/paf)
  <ASM2>  asm2 alignment file (sam/bam/cram/paf)

Options:
  -1, --s1 <NAME>           label for asm1 sample (used in output names) [default: asm1]
  -2, --s2 <NAME>           label for asm2 sample (used in output names) [default: asm2]
      --paf                 input files are PAF
      --ms                  use ms:i: tag rather than AS:i: for alignment score
  -p, --partition           write one file per haplotype instead of a single merged output file
  -b, --both                write reads with equal alignment scores to both output files (requires -p)
  -k, --keep-loser          also write each read's best alignment to the losing haplotype as a secondary record, tagged hs:A: (merged output only)
      --loser-frac <FLOAT>  score floor for AS score fraction in --keep-loser [default: 0.8]
  -o, --output <FILE>       output file name [default: hiphap_{s1}_{s2}_merged.*]
  -u, --unmapped <DEST>     where to write reads unmapped in both assemblies: asm1, asm2, or discard [default: asm1] [possible values: asm1, asm2, discard]
      --ref-merged <FILE>   combined reference FASTA for merged CRAM output; required for merged CRAM output
      --ref1 <FILE>         reference FASTA for cram file (asm1)
      --ref2 <FILE>         reference FASTA for cram file (asm2)
  -A, --match-sc <FLOAT>    per-base match score from aligner scoring scheme (auto-estimated if omitted)
      --no-hapq             skip HAPQ score calculation and hq tag output (e.g. for comparing GRCh38 vs CHM13)
      --no-span-chrom       disable writing the chromosome-spanning reads file (*_span_chrom.fastq, or .txt for PAF)
  -t, --threads <INT>       number of threads[default: 6; 8 with -p]
  -h, --help                Print help
  -V, --version             Print version
```

Each output record is annotated with two tags: an `HP:i:` tag carrying the haplotype it was assigned to (`1` = asm1, `2` = asm2), and an `hq:i:` tag carrying the HapQ score, unless `--no-hapq` is set. See [Haplotype tag (HP)](#haplotype-tag-hp).

Pass `--keep-loser` to also write each read's alignment to the *losing* haplotype, flagged secondary and marked with an `hs:A:` tag — see [Losing-haplotype alignments](#losing-haplotype-alignments---keep-loser).

The per-base match score used by the HAPQ calculation is auto-estimated from the `ms:i:` tags of the input files, and the estimate is printed at the start of the run. Pass `-A`/`--match-sc <FLOAT>` to set it explicitly — it must match the aligner's `-A`, which differs by preset: `1` for minimap2 `map-hifi`, `2` for `map-ont` and `map-pb`. Auto-estimation fails if the inputs carry no usable `ms:i:` tags, in which case `-A` is required.

### Threads

`-t`/`--threads` is a total budget, divided between the two input readers (one per assembly) and the output writers (one merged writer, or one per haplotype with `-p`). For BAM and CRAM a writer is weighted more heavily than a reader, since compressing a record costs several times as much as decompressing one: **4x** for a merged run, where a single writer handles both assemblies' output, and **3x** with `-p`, where the work is split across two writers. For uncompressed SAM the weight is 1x. That makes the smallest fully balanced budget 6 for a merged run (`1 + 1` readers, `4` writer) and 8 with `-p` (`1 + 1` readers, `3 + 3` writers), which are the two defaults.

Larger budgets keep the same ratio, so 10 merged gives `2 + 2` readers and `6` writer, and 16 with `-p` gives `2 + 2` readers and `6 + 6` writers. A budget too small to give every file one thread is raised to the minimum, and in `-p` mode an odd leftover thread is left unused so the two haplotype writers stay symmetric. Note that htslib spawns these as *extra* worker threads alongside the main thread, so the process peaks at roughly `-t` plus one.

`--threads` has no effect in PAF mode, which is single-threaded text I/O; passing it there prints a warning. `--ref1`/`--ref2`/`--ref-merged` are likewise ignored with `--paf`.

## Example Workflow

HipHap reads the two inputs as parallel streams and pairs records by read, so **both files must present reads in the same order**. Raw [minimap2](https://github.com/lh3/minimap2) output satisfies this (`SO:unsorted GO:query` — records grouped by read, in input order), as does name-sorting both files. A coordinate-sorted input does not, and neither does name-sorting one file while leaving the other in aligner order: the run aborts with `alignment streams out of sync` rather than producing wrong output.

```bash
samtools sort -n -o name_sort.bam coord_sort.bam   # apply to BOTH inputs
```

Truncating an input to a record count (`head -n`) also breaks this, since the cut lands mid-read and at a different read in each file; subset by read name instead.

### Diploid assembly alignment

```bash
# If needed, split diploid genome assembly FASTA into respective haplotypes
separate_haps_fasta -1 MATERNAL -2 PATERNAL hg002v1.1.fa
# writes hg002v1.1.MATERNAL.fa and hg002v1.1.PATERNAL.fa

# Align reads to each haplotype
minimap2 -ax map-hifi -o asm1_alignments.sam hg002v1.1.MATERNAL.fa reads.fastq
minimap2 -ax map-hifi -o asm2_alignments.sam hg002v1.1.PATERNAL.fa reads.fastq

# Run hiphap
hiphap -1 mat -2 pat asm1_alignments.sam asm2_alignments.sam
# Output: hiphap_mat_pat_merged.sam  hiphap_mat_pat_span_chrom.fastq

# Save as sorted BAM
samtools sort -@ 12 -o hiphap_mat_pat_merged.bam hiphap_mat_pat_merged.sam
```

By default HipHap writes a **single merged output file** with a merged header, so no separate `samtools merge` step is needed. Each record carries an `HP:i:1`/`HP:i:2` tag naming the haplotype it was assigned to (see [Haplotype tag (HP)](#haplotype-tag-hp)). Use `-o` to name it:

```bash
hiphap -1 mat -2 pat -o sample.sam asm1_alignments.sam asm2_alignments.sam
# Output: sample.sam  sample_span_chrom.fastq
```

HipHap also writes `hiphap_{s1}_{s2}_span_chrom.fastq`, a FASTQ of reads whose alignments span more than one chromosome. These reads are emitted for easy realignment. Pass `--no-span-chrom` to skip writing this file. If inputs are PAF files, a tsv of the chromosome spamnign reads rather than a fastq is output, as the read sequences are not stored in the input PAF files. When `-o` is given, this file is named after it (`sample_span_chrom.fastq` above) so it lands beside the alignment output.

#### Notes:
- Merged output requires the **two assemblies to have unique contig names**. The merged header concatenates the two inputs `@SQ` lists, so any contig name shared between them is ambiguous; pass `-p` for such inputs.
- The output format always follows the input format. An `-o` extension that disagrees is ignored, with a warning.
- For **CRAM** inputs, merged output is also CRAM and requires a single combined reference (`--ref-merged <FILE>`) containing all contigs of both haplotypes, typically the original diploid assembly:

  ```bash
  hiphap --ref1 hg002v1.1.MATERNAL.fa --ref2 hg002v1.1.PATERNAL.fa \
         --ref-merged hg002v1.1.fa mat_alignments.cram pat_alignments.cram
  ```

### One file per haplotype (`-p`)

`-p`/`--partition` writes one output file per haplotype instead of the merged file:

```bash
hiphap -p -1 mat -2 pat asm1_alignments.sam asm2_alignments.sam
# Output: hiphap_mat.sam  hiphap_pat.sam  hiphap_mat_pat_span_chrom.fastq

# Merge them afterwards if desired
samtools merge -@ 12 merged.sam hiphap_mat.sam hiphap_pat.sam
```

With `-p`, `-o` is the stem the two files share:

```bash
hiphap -p -1 mat -2 pat -o sample asm1_alignments.sam asm2_alignments.sam
# Output: sample_mat.sam  sample_pat.sam  sample_span_chrom.fastq
```

#### Notes:
- `--both` requires `-p`: in a merged file a tied read would get two primary records.
- `--keep-loser` cannot be used with `-p`: a per-haplotype file holds one assembly's contigs, so there is nowhere in it to put the other haplotype's alignments. See [Losing-haplotype alignments (`--keep-loser`)](#losing-haplotype-alignments---keep-loser).
- Contig names may be shared between the inputs in this mode, since each haplotype keeps its own header.

### Comparing different reference genomes

HipHap can also be used to select best alignment of a read between different reference genomes (e.g. GRCh38 and CHM13). For this use case, the HAPQ score is generally not meaningful, so pass `--no-hapq` to skip its calculation:

```bash
hiphap --no-hapq -p -1 grch38 -2 chm13 grch38_alignments.sam chm13_alignments.sam
# Output: hiphap_grch38.sam  hiphap_chm13.sam
```

`-p` is needed here because the two references share contig names (`chr1` in both), which merged output cannot represent.

### CRAM input files

If input files are CRAM format, the original reference genomes must be provided. Merged output additionally needs a combined reference (`--ref-merged`), since a single CRAM writer needs one reference covering all contigs of both haplotypes:

```bash
hiphap --ref1 asm1_hap.fasta --ref2 asm2_hap.fasta --ref-merged diploid.fasta \
       asm1_alignments.cram asm2_alignments.cram

# per-haplotype output needs only the two originals
hiphap -p --ref1 asm1_hap.fasta --ref2 asm2_hap.fasta asm1_alignments.cram asm2_alignments.cram
```

## Example PAF Usage

 #### Notes: 
- It is important to use the `--paf-no-hit` flags when aligning with minimap2 to output unmapped reads to the file
- If a SAM file is converted to a PAF file with `paftools.js sam2paf`, it will **NOT** have the required AS:i: tag and HipHap will fail to run

```bash
minimap2 -cx map-hifi --paf-no-hit -o asm1_alignments.paf hg002v1.1.MATERNAL.fa reads.fastq
minimap2 -cx map-hifi --paf-no-hit -o asm2_alignments.paf hg002v1.1.PATERNAL.fa reads.fastq

hiphap --paf out_asm1.paf out_asm2.paf
# Output: hiphap_asm1_asm2_merged.paf  hiphap_asm1_asm2_span_chrom.txt
```

`-p` also works with PAF

```bash
hiphap --paf -p out_asm1.paf out_asm2.paf
# Output: hiphap_asm1.paf  hiphap_asm2.paf  hiphap_asm1_asm2_span_chrom.txt
```

## Weighted Alignment Scoring Mechanism

For each read, HipHap computes a single weighted alignment score per assembly using all primary and supplementary alignments (secondary alignments are passed through to the output but ignored when scoring).

Let $n$ be the number of primary and supplementary alignments for a given read to that reference genome.

Let $L$ be the full read length (sum of query-consuming CIGAR operations on a non-secondary record, or the `qlen` field of the PAF record).

Let $I_i = [r_i^{\mathrm{start}}, r_i^{\mathrm{end}})$ denote the interval on the read covered by alignment $i$.

Let

$$
B = \left| \bigcup_{i=1}^{n} I_i \right|
$$

be the number of **read bases** covered by at least one alignment.

Let $a_i$ be the alignment score for alignment $i$ (`AS:i:` by default, or `ms:i:` if `--ms` is set).

Let $l_i$ be the alignment length in read coordinates for alignment $i$.

$$
S = \frac{\sum_{i=1}^{n} a_i}{\sum_{i=1}^{n} l_i} \cdot B \cdot \frac{B}{L}
$$

The first factor is the average alignment score per aligned base. It is multiplied by the number of unique read bases covered, $B$, and then scaled by the read coverage fraction $B/L$, so that reads which align over a large fraction of their length are weighted more heavily than reads which align only over a small portion.

For each read, the assembly with the higher $S$ wins; its full alignment cluster (including secondary alignments) is written to the corresponding output file. If $S$ is equal in both assemblies, the "better" assignment is determined by a hash of the read name, or the read is written to both output files when `--both` is used.

The ratio of the two $S$ values is also what `--loser-frac` thresholds on when [keeping the losing haplotype's alignments](#losing-haplotype-alignments---keep-loser).

## Haplotype tag (HP)

Every assigned read carries an `HP:i:` tag recording **which haplotype HipHap chose**:

| Tag | Meaning |
|---|---|
| `HP:i:1` | assigned to asm1 (`-1`/`--s1`) |
| `HP:i:2` | assigned to asm2 (`-2`/`--s2`) |
| *absent* | read unmapped in both assemblies (no assignment was made) |

This is the same tag convention used by phasing tools such as WhatsHap and HapCUT2, so IGV groups and colors reads by haplotype out of the box, and the assignment can be recovered directly rather than by parsing contig names:

```bash
samtools view -d HP:1 hiphap_mat_pat_merged.bam   # only reads assigned to asm1
```

In merged output the two `@SQ` lists are concatenated, so without this tag the assignment is only recoverable from the contig a read landed on, which requires the assembly to use a recognizable per-haplotype contig naming convention. `HP` is written in every mode (merged, `-p`, and PAF) and is independent of `--no-hapq`.

For reads tied between the two assemblies the `HP` tag reflects whichever side the read was written to; such reads are identified by `hq:i:0` (and appear in both outputs, with `HP:i:1` and `HP:i:2` respectively, under `-p --both`).

## Losing-haplotype alignments (`--keep-loser`)

By default only the winning haplotype's alignments are written and the losing haplotype's are discarded. When the two haplotypes scored close to each other, that discarded alignment is the evidence for *why* the call was marginal, and it is what you would want to see in IGV next to the winner. `--keep-loser` also writes it:

```bash
hiphap --keep-loser -1 mat -2 pat asm1_alignments.bam asm2_alignments.bam
```

For each read, the losing haplotype's **primary and supplementary** alignments are appended to the merged file directly after the winning cluster, rewritten as secondary alignments. The losing haplotype's own secondary alignments are dropped, as is an unmapped record, since neither says anything about the call.

Each such record carries an `hs:A:` tag — *haplotype status* — naming what the record was in its own cluster, which is no longer recoverable from the rewritten `FLAG`:

| Tag | Meaning | Resulting FLAG |
|---|---|---|
| `hs:A:P` | was the **primary** alignment to the losing haplotype | `0x100` (secondary) |
| `hs:A:S` | was a **supplementary** alignment to the losing haplotype | `0x900` (secondary + supplementary) |
| *absent* | **winner** — this record belongs to the haplotype the read was assigned to | unchanged |

`hs` is the only field that identifies a losing record. Every read has exactly one `hs:A:P`, so that is also what to count:

```bash
samtools view -d hs   merged.bam    # every losing-haplotype record
samtools view -d hs:P merged.bam    # best alignment to the losing haplotype, one per read
samtools view -d hs:S merged.bam    # the losing cluster's split segments
samtools view -e '![hs]' merged.bam # winners only — the file hiphap writes without the flag
```

Two things that look like they would identify a losing record but do **not**:

- **`FLAG` `0x100`** is also carried by the winning haplotype's own secondary alignments, which are passed through unchanged. Selecting on `-f 0x100` mixes the two. Filtering with `-F 0x900` is not "winners only" either — it additionally drops the winner's supplementary and secondary records, so it returns only winning *primary* alignments rather than everything hiphap writes without `--keep-loser`.
- **`tp:A:`** is minimap2's tag, written before hiphap ran, and records whether an alignment was primary *in its own input file*. Every losing record was primary in its own file, so losers keep `tp:A:P` and cannot be separated from winners by it.

Setting `0x100` rather than leaving these records primary keeps the invariant that each read has exactly one primary alignment in the merged file, so read counts, coverage, and consensus tools are unaffected by the flag. The `0x800` bit is left in place, so a former supplementary is `0x900` and the split-alignment structure survives.

Details:

- **`SEQ` and `QUAL` are `*`.** The winning record for the same read already holds the read's bases, and writing them a second time would roughly double the file. This is what minimap2 itself does for secondary alignments.
- **`hq:i:` is written** with the same HapQ as the winning record, so the confidence in the call is visible on both sides of it (omitted under `--no-hapq`, as elsewhere).
- **No `HP:i:` tag.** `HP` records the haplotype hiphap *assigned* the read to, and these records are the haplotype it was assigned away from. Tagging them would put them in the wrong bin for any tool that splits on `HP`.
- **`MAPQ` is unchanged** — it is still the mapping quality of that alignment within the losing haplotype.
- These records do not appear in the chromosome-spanning reads file, which tracks winning clusters only.
- `--keep-loser` cannot be combined with `-p`/`--partition`.

### Score threshold (`--loser-frac`)

A clear loser is not informative, so a losing alignment is only written when its weighted score (the score that picked the winner, see [Weighted Alignment Scoring Mechanism](#weighted-alignment-scoring-mechanism)) is at least `--loser-frac` of the winner's — **0.8 by default**. Reads tied between the haplotypes (`hq:i:0`) always qualify.

Raise it to keep only the closest calls, or lower it to keep more:

```bash
hiphap --keep-loser --loser-frac 0.95 asm1.bam asm2.bam   # near-ties only
hiphap --keep-loser --loser-frac 1.0  asm1.bam asm2.bam   # exact ties only
```

The value must be in `(0, 1]`; `0` is rejected rather than treated as "keep everything", since a loser scoring nothing is not evidence. At `1.0` only exact ties survive, because a clear winner by definition scores above its loser.

The run summary reports how many reads cleared the threshold:

```
Reads with losing-hap alignments kept (>= 0.80 of winner): 17897 (91.9% of assigned reads)
```

### PAF

`--keep-loser` works for merged PAF output too. Losing lines get their `tp:A:P` rewritten to `tp:A:S` and the `hq:i:`/`hs:A:` tags appended, and lines already marked `tp:A:S` or with a `*` target are dropped.

In PAF `hs` is always `P`: minimap2 gives every non-secondary chain `tp:A:P`, including the ones that become supplementary in SAM, so a PAF line carries nothing that separates the two. `hs:A:S` only ever appears in SAM/BAM/CRAM output.

## Haplotype Assignment Quality (HapQ)

For each read assigned to a winning haplotype, HipHap reports a HapQ score in the `hq:i:` tag of the output record. HapQ is a Phred-like confidence [0-60] that the read was assigned to the correct haplotype. The calculation is modeled on BWA-MEM's [`mem_approx_mapq_se`](https://github.com/lh3/bwa/blob/b92993c1161e73167181558856567ef2f367e3f0/bwamem.c#L982-L1006).

Let $S_w$ and $S_l$ be the weighted alignment scores of the winning and losing assemblies, $m$ the per-base match score of the alignment software used (`-A`/`--match-sc`), and $k$ the number of non-secondary alignments (splits) on the winning side.

HapQ is the product of: 

d, the approximate difference, in matching bases, between the winning and losing alignments between haplotypes.
```math
d = \frac{S_w - S_l}{m}
```

The split penalty, $\rho$ , which penalizes reads with more than three split alignments, which often fall in complex/repetitive regions where haplotype assignment is less reliable.
```math
\rho = \begin{cases} 1 & k \le 3 \\ \frac{3}{k} & k > 3 \end{cases}
```

The final HapQ calculation is then: 
```math
\text{HapQ}\;=\;\lfloor\, 6.02 \cdot \rho \cdot d \,\rfloor
```

Clamped to the range [0,60]. 

Special cases:
- Read mapped in only one assembly: HapQ = 60.
- Read tied between assemblies (winner = `Both`): HapQ = 0.
- Read unmapped in both assemblies: no `hq` tag is written.
- A non-tied read whose score rounds below 1 is reported as HapQ = 1, not 0, so that `hq:i:0` means a true tie and nothing else. HapQ = 0 is therefore never the result of rounding.

If `--no-hapq` is set, HAPQ is not computed and no `hq:i:` tag is added (recommended when the two inputs are not haplotypes of the same sample, e.g. GRCh38 vs CHM13).

## Citation
If HipHap has helped you in your research, please cite our preprint at: TODO
