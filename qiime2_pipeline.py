#!/usr/bin/env python3
"""
Metagenomic Analysis using QIIME2
Aim: Assess microbial community structure, diversity, and taxonomic
     composition from 16S amplicon sequencing data.

Requires: QIIME 2 installed (e.g. in a conda environment) and the `qiime`
          command available on $PATH, plus a pre-trained taxonomy classifier
          (e.g. gg-13-8-99-515-806-nb-classifier.qza)
"""

import subprocess

METADATA = "sample-metadata.tsv"
CLASSIFIER = "gg-13-8-99-515-806-nb-classifier.qza"
SAMPLING_DEPTH = "1103"   # chosen from the feature-table summary (table.qzv)


def run(cmd):
    print(f">>> {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def import_sequences():
    run([
        "qiime", "tools", "import",
        "--type", "EMPSingleEndSequences",
        "--input-path", "input_qiime",
        "--output-path", "single-end-sequences.qza",
    ])


def demultiplex():
    run([
        "qiime", "demux", "emp-single",
        "--i-seqs", "single-end-sequences.qza",
        "--m-barcodes-file", METADATA,
        "--m-barcodes-column", "barcode-sequence",
        "--o-per-sample-sequences", "demux.qza",
        "--o-error-correction-details", "demux-details.qza",
    ])
    run([
        "qiime", "demux", "summarize",
        "--i-data", "demux.qza",
        "--o-visualization", "demux.qzv",
    ])


def denoise_dada2():
    run([
        "qiime", "dada2", "denoise-single",
        "--i-demultiplexed-seqs", "demux.qza",
        "--p-trim-left", "0",
        "--p-trunc-len", "120",
        "--o-representative-sequences", "rep-seqs.qza",
        "--o-table", "table.qza",
        "--o-denoising-stats", "stats.qza",
    ])
    run([
        "qiime", "metadata", "tabulate",
        "--m-input-file", "stats.qza",
        "--o-visualization", "stats.qzv",
    ])


def summarize_feature_table():
    run([
        "qiime", "feature-table", "summarize",
        "--i-table", "table.qza",
        "--m-sample-metadata-file", METADATA,
        "--o-visualization", "table.qzv",
    ])
    run([
        "qiime", "feature-table", "tabulate-seqs",
        "--i-data", "rep-seqs.qza",
        "--o-visualization", "rep-seqs.qzv",
    ])


def build_phylogenetic_tree():
    run([
        "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
        "--i-sequences", "rep-seqs.qza",
        "--o-alignment", "aligned-rep-seqs.qza",
        "--o-masked-alignment", "masked-aligned-rep-seqs.qza",
        "--o-tree", "unrooted-tree.qza",
        "--o-rooted-tree", "rooted-tree.qza",
    ])
    run(["qiime", "tools", "export", "--input-path", "rooted-tree.qza", "--output-path", "rooted-tree"])


def diversity_analysis():
    run([
        "qiime", "diversity", "core-metrics-phylogenetic",
        "--i-phylogeny", "rooted-tree.qza",
        "--i-table", "table.qza",
        "--p-sampling-depth", SAMPLING_DEPTH,
        "--m-metadata-file", METADATA,
        "--output-dir", "diversity-core-metrics-phylogenetic",
    ])

    # Alpha diversity significance (Faith's PD, evenness)
    run([
        "qiime", "diversity", "alpha-group-significance",
        "--i-alpha-diversity", "diversity-core-metrics-phylogenetic/faith_pd_vector.qza",
        "--m-metadata-file", METADATA,
        "--o-visualization", "faith-pd-group-significance.qzv",
    ])
    run([
        "qiime", "diversity", "alpha-group-significance",
        "--i-alpha-diversity", "diversity-core-metrics-phylogenetic/evenness_vector.qza",
        "--m-metadata-file", METADATA,
        "--o-visualization", "evenness-group-significance.qzv",
    ])

    # Beta diversity significance (unweighted UniFrac, by body-site / subject)
    for column, out_name in [
        ("body-site", "unweighted-unifrac-body-site-group-significance.qzv"),
        ("subject", "unweighted-unifrac-subject-group-significance.qzv"),
    ]:
        run([
            "qiime", "diversity", "beta-group-significance",
            "--i-distance-matrix",
            "diversity-core-metrics-phylogenetic/unweighted_unifrac_distance_matrix.qza",
            "--m-metadata-file", METADATA,
            "--m-metadata-column", column,
            "--p-pairwise",
            "--o-visualization", out_name,
        ])


def emperor_plot():
    run([
        "qiime", "emperor", "plot",
        "--i-pcoa", "diversity-core-metrics-phylogenetic/unweighted_unifrac_pcoa_results.qza",
        "--m-metadata-file", METADATA,
        "--p-custom-axes", "days-since-experiment-start",
        "--o-visualization", "unweighted-unifrac-emperor-days-since-experiment-start.qzv",
    ])


def alpha_rarefaction():
    run([
        "qiime", "diversity", "alpha-rarefaction",
        "--i-table", "table.qza",
        "--i-phylogeny", "rooted-tree.qza",
        "--p-max-depth", "4000",
        "--m-metadata-file", METADATA,
        "--o-visualization", "alpha-rarefaction.qzv",
    ])


def taxonomic_classification():
    run([
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-classifier", CLASSIFIER,
        "--i-reads", "rep-seqs.qza",
        "--o-classification", "taxonomy.qza",
    ])
    run([
        "qiime", "metadata", "tabulate",
        "--m-input-file", "taxonomy.qza",
        "--o-visualization", "taxonomy.qzv",
    ])
    run([
        "qiime", "taxa", "barplot",
        "--i-table", "table.qza",
        "--i-taxonomy", "taxonomy.qza",
        "--m-metadata-file", METADATA,
        "--o-visualization", "taxa-bar-plots.qzv",
    ])


def main():
    import_sequences()
    demultiplex()
    denoise_dada2()
    summarize_feature_table()
    build_phylogenetic_tree()
    diversity_analysis()
    emperor_plot()
    alpha_rarefaction()
    taxonomic_classification()

    print("QIIME2 pipeline complete. Key outputs: table.qzv, taxa-bar-plots.qzv,")
    print("  alpha-rarefaction.qzv, unweighted-unifrac-*-group-significance.qzv")


if __name__ == "__main__":
    main()
