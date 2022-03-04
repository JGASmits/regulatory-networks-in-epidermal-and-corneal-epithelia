import os
import pandas as pd
import subprocess as sp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import csv
from pathlib import Path
import qnorm    
import math

def distance_weight(include_promoter=False, include_enhancer=True, alpha=1e5, padding=100000, keep1=5000, remove=2000):
    """
    Built weight distribution from TSS. (From Ananse)
    """
    # alpha is half site, default setting is 1e5, which means at 1e5 position weight is 0.5
    # padding is the full range we used
    # remove is promoter removed range
    # keep1 is keep full binding score range

    u = -math.log(1.0 / 3.0) * 1e5 / alpha

    if include_promoter and include_enhancer:
        weight1 = pd.DataFrame(
            {"weight": [1 for z in range(1, remove + 1)], "dist": range(1, remove + 1)}
        )
        weight2 = pd.DataFrame(
            {
                "weight": [1 for z in range(remove + 1, keep1 + 1)],
                "dist": range(remove + 1, keep1 + 1),
            }
        )
        weight3 = pd.DataFrame(
            {
                "weight": [
                    2.0
                    * math.exp(-u * math.fabs(z) / 1e5)
                    / (1.0 + math.exp(-u * math.fabs(z) / 1e5))
                    for z in range(1, padding - keep1 + 1)
                ],
                "dist": range(keep1 + 1, padding + 1),
            }
        )
    elif not include_promoter and include_enhancer:
        weight1 = pd.DataFrame(
            {"weight": [0 for z in range(1, remove + 1)], "dist": range(1, remove + 1)}
        )
        weight2 = pd.DataFrame(
            {
                "weight": [1 for z in range(remove + 1, keep1 + 1)],
                "dist": range(remove + 1, keep1 + 1),
            }
        )
        weight3 = pd.DataFrame(
            {
                "weight": [
                    2.0
                    * math.exp(-u * math.fabs(z) / 1e5)
                    / (1.0 + math.exp(-u * math.fabs(z) / 1e5))
                    for z in range(1, padding - keep1 + 1)
                ],
                "dist": range(keep1 + 1, padding + 1),
            }
        )
    elif include_promoter and not include_enhancer:
        weight1 = pd.DataFrame(
            {"weight": [1 for z in range(1, remove + 1)], "dist": range(1, remove + 1)}
        )
        weight2 = pd.DataFrame(
            {
                "weight": [0 for z in range(remove + 1, keep1 + 1)],
                "dist": range(remove + 1, keep1 + 1),
            }
        )
        weight3 = pd.DataFrame(
            {
                "weight": [0 for z in range(1, padding - keep1 + 1)],
                "dist": range(keep1 + 1, padding + 1),
            }
        )
    else:
        weight1 = pd.DataFrame(
            {"weight": [0 for z in range(1, remove + 1)], "dist": range(1, remove + 1)}
        )
        weight2 = pd.DataFrame(
            {
                "weight": [0 for z in range(remove + 1, keep1 + 1)],
                "dist": range(remove + 1, keep1 + 1),
            }
        )
        weight3 = pd.DataFrame(
            {
                "weight": [0 for z in range(1, padding - keep1 + 1)],
                "dist": range(keep1 + 1, padding + 1),
            }
        )

    weight = pd.concat([weight1, weight2, weight3])
    #weight_dict = pd.Series(weight.weight.values, index=weight.dist).to_dict()
    return weight



def TSS_window_to_region(genome_TSS_gtf, bedregion_file, TSS_window):
    """
    Takes a genome TSS gtf file or bedfile as input (see the 'genome_TSS_GTF' script) together with a bed file with regions and
    a window length. 
    It wil subsequently take a window around the TSS, and list all regions overlapping with this window. It will
    furhtermore calculate the distance of te bedregion to the TSS. 
    It will return a dataframe containing all these results.
    """
    sp.check_call(
        f"nice -14 bedtools window "
        f"-w {TSS_window} "
        f"-a {genome_path_gtf_small} "
        f"-b {bedregion_file}  "
        f"> tmp_bedtool_window.txt",
        shell=True,
    )

    pd_table = pd.read_table(
        "tmp_bedtool_window.txt",
        index_col=None,
        names=[
            "chr_TSS",
            "start_TSS",
            "end_TSS",
            "gene_name",
            "Chrom",
            "ChromStart",
            "ChromEnd",
            "region_annotation",
        ],
        dtype={
            "chr_TSS": str,
            "start_TSS": int,
            "end_TSS": int,
            "gene_name": str,
            "Chrom": str,
            "ChromStart": int,
            "ChromEnd": int,
            "region_annotation": str,
        },
    )

    os.remove("tmp_bedtool_window.txt")

    distance = (
        []
    )  # calculate the distance based on the TSS start end region start site, check which one is upstream for distance calc.
    for x, y in zip(pd_table["start_TSS"], pd_table["ChromStart"]):
        if x > y:
            distance.append(x - y)
        elif y > x:
            distance.append(y - x)
        else:
            distance.append("NaN")

    pd_table["distance"] = distance
    return pd_table

def make_bedfile_from_column(df, column_name, filename):
    """
    df: a df containing a column with : chr1:100-1000 (or 1:100-1000)
    column_name: the name of the column with the loc ids that need to be exported to a bed file
    filename: bedfile the location will be exported to.
    """
    df[["chrom", "chromStart", "chromEnd"]] = df.loc[:,column_name].str.split(
       "[:-]", expand=True
     )
    df = df.loc[:, ["chrom", "chromStart", "chromEnd", 'annotation']]
    df.to_csv(
    filename, sep="\t", header=False, index=False
)

def make_chrom_chromStart_chromEnd_bedfile(df, filename, column_name="index"):
    """
    df: a df containing a column with : chr1:100-1000 (or 1:100-1000)
    column_name: either index or the name of the column with the loc ids that need to be exported to a bed file
    filename: bedfile the location will be exported to.
    """
    if column_name == "index":
        df[["chrom", "chromStart", "chromEnd"]] = df.index.str.split(
            "[:-]", expand=True
        )
    else:
        df[["chrom", "chromStart", "chromEnd"]] = df.loc[:, column_name].str.split(
            "[:-]", expand=True
        )

    df = df.loc[:, ["chrom", "chromStart", "chromEnd"]]
    df.to_csv(filename, sep="\t", header=False, index=False)

def TSS_to_region(genome_TSS_gtf, bedregion_file, extra_options):
    """s
    Takes a genome TSS gtf file or bedfile as input (see the 'genome_TSS_GTF' script) together with a bed file with regions. It wil
    subsequently map the regions from the bed file to the closest TSS region
    The extra_options is a aditional string with options for bedtool closest: e.g. '-k 1 -d ' 
    """
    sp.check_call(
        f"nice -14 bedtools closest "
        f"{extra_options} "
        f"-a {genome_TSS_gtf} "
        f"-b {bedregion_file}  "
        f"> tmp_bedtool_closest.txt",
        shell=True,
    )

    pd_table = pd.read_table(
        "tmp_bedtool_closest.txt",
        index_col=None,
        names=[
            "chr_TSS",
            "start_TSS",
            "end_TSS",
            "gene_name",
            "Chrom",
            "ChromStart",
            "ChromEnd",
            "gene_region_dist",
        ],
        dtype={
            "chr_TSS": str,
            "start_TSS": int,
            "end_TSS": int,
            "gene_name": str,
            "Chrom": str,
            "ChromStart": int,
            "ChromEnd": int,
            "region_annotation": str,
            "gene_region_dist": int,
        },
    )

    os.remove("tmp_bedtool_closest.txt")

    return pd_table

def select_DEGS_from_DESEQ_results(
    deseq_file, pval_cutoff=0.05, log2FC_cutoff=0.68, min_counts=100
):
    """
    select all DEGS from a deseq2 result csv file, with a aditional column called total_counts
    that contains the total counts of each gene.
    Use this object to annotate all genes that pass the pval_cutoff, FC_cutoff and minimum ammount
    of counts needed.
    Add a column called 'gene annotation' with a 1 for
    all genes that are differential, and a 0 for all nondiff genes
    """

    DEG_pd = pd.read_csv(
        deseq_file,
        sep=",",
        index_col=0,
        dtype={
            "baseMean": float,
            "log2FoldChange": float,
            "lfcSE": float,
            "stat": float,
            "pvalue": float,
            "padj": float,
            "gene": str,
            "total_counts": float,
        },
    )

    DEG_pd["gene_annotation"] = (
        (DEG_pd.log2FoldChange.astype(float).abs() >= log2FC_cutoff)
        & (DEG_pd.padj.astype(float).abs() <= pval_cutoff)
        & (DEG_pd.total_counts.astype(int) > min_counts)
    )

    DEG_pd["gene_annotation"] = DEG_pd["gene_annotation"].astype(int)

    return DEG_pd


def bedtool_closest_annotated(genome_TSS_gtf, bedregion_file):
    """
    Takes a genome TSS gtf file as input (see the 'genome_TSS_GTF' script) together with a bed file with a 4th annotation column and
    with regions. It wil
    subsequently map the regions from the bed file to the closest TSS region
    """
    sp.check_call(
        f"nice -14 bedtools closest "
        f"-k 1 "  # report the clossest gene
        f"-d "  # reports the distance to the gene
        f"-b {genome_TSS_gtf} "
        f"-a {bedregion_file}  "
        f"> tmp_bedtool_closest.txt",
        shell=True,
    )

    pd_table = pd.read_table(
        "tmp_bedtool_closest.txt",
        index_col=None,
        names=[
            "chr",
            "start",
            "end",
            "annotation",
            "gene_chr",
            "gene_start",
            "gene_end",
            "gene",
            "region-gene-dist",
        ],
        dtype={
            "chr_e": str,
            "start_e": int,
            "end_e": int,
            "chr": str,
            "start": int,
            "end": int,
            "annotation": str,
            "gene": str,
            "enh-gene-dist": str,
        },
    )

    os.remove("tmp_bedtool_closest.txt")

    return pd_table

def hello_world(chr_list1, db_location1):
    """
    chr_list = "list of chr identifiers, or chr locations (either "1" or "1:1002-3999") that need
    to be converted (into "chr1" or "chr1:1002-3999")
    db_location = "tab seperated identifier linking file"
    """
    table = pd.read_table(db_location1)
    table = dict(zip(table.ensembl, table.UCSC))

    new_ids = []

    if all(chr_list1.str.contains(":")):
        print("list contains the ':' symbol, replacing chrom ids in front of :")
        for line in chr_list1:
            chrom = line.split(":")[0]
            coordinate = line.split(":")[1]
            chrom = table.get(chrom)
            new_ids.append(chrom + ":" + coordinate)
    else:
        print("no ':' symbol found, asuming only 1 chrom id per line")
        for line in chr_list1:
            new_ids.append(table.get(line))

    return new_ids

def summits_2_regions(summit_file, output_regions_file, window_size):
    """
    generates a window across the summits of a bed file
    summit_file = "bed file containing summit regions"
    output_regions_file = 'file location the extended summits need to be exported to'
    window_size = 'total length of the window spanning the summit'
    """
    summits = bed2df(summit_file)
    summits["enhChrom"] = summits["chrom"]
    summits["enhStart"] = (summits.chromStart - window_size / 2).astype(int)
    summits["enhEnd"] = (summits.chromStart + window_size / 2).astype(int)
    summits.loc[:, ["enhChrom", "enhStart", "enhEnd"]].to_csv(
        output_regions_file, sep="\t", header=False, index=False
    )

def coverage_table_normalization(dataframe):
    """
    Normalizes a df from a coverage table, first log10 transforming followed by quantile normalization.
    dataframe = "input dataframe containing the output of the covarage table function"
    """
    # log10 transform the reads in the coverage table
    log_df = np.log10(dataframe + 1)
    # perform Quantile normalization:
    qnorm_df = qnorm.quantile_normalize(log_df)
    #rank_mean = (
    #    log_df.stack().groupby(log_df.rank(method="first").stack().astype(int)).mean()
    #)

    #norm_df = log_df.rank(method="min").stack().astype(int).map(rank_mean).unstack()
    #return norm_df
    return(qnorm_df)

def make_bigbed_UCSC_genome_browser(
    loc_list, output_file, output_bigbed_file, db_location, genome_path_size_UCSC
):
    """
    chr_list = "list of chr identifiers, or chr locations (either "1" or "1:1002-3999") that need 
    to be visualized on the UCSC genome browser
    output_file = 'location to save the regular bedfile to'
    output_bigbed_file = "location to save the bigbed file to"
    db_location = 'location of a table containing UCSC 2 ensamble chrom ID's, if using Ensample mapped data
    'this is needed to change the chrom IDs'
    """
    chrom_table = pd.DataFrame(loc_list)
    # make a dictionary containing Ensembl ids as keys and UCSC ids as vallues
    table = pd.read_table(db_location)
    table = dict(zip(table.ensembl, table.UCSC))

    # load the chrom Ids
    chrom_table[["chrom", "start", "end"]] = chrom_table["loc"].str.split(
        "[:-]", expand=True
    )
    chrom_table = chrom_table.astype({"chrom": str, "start": int, "end": int})

    # generate the bed file output (ensemble IDs)
    np.savetxt(
        output_file, chrom_table[["chrom", "start", "end"]], delimiter="\t", fmt="%s"
    )
    sp.check_call(f"bedSort {output_file} {output_file}", shell=True)

    # map the chrom id's  to UCSC ids (necisary for UCSC browser vizualization)
    chrom_table["chrom"] = chrom_table["chrom"].map(table)
    chrom_table["loc"] = (
        chrom_table["chrom"].astype(str)
        + ":"
        + chrom_table["start"].astype(str)
        + "-"
        + chrom_table["end"].astype(str)
    )
    chrom_table = chrom_table.drop(["loc"], axis=1)

    # save the chrom loc to a bed file:
    np.savetxt("./tmp.txt", chrom_table, delimiter="\t", fmt="%s")
    sp.check_call(f"bedSort ./tmp.txt ./tmp.txt", shell=True)
    sp.check_call(
        f"bedToBigBed -type=bed3 ./tmp.txt {genome_path_size_UCSC} {output_bigbed_file}",
        shell=True,
    )
    sp.check_call(f"rm -r ./tmp.txt", shell=True)

def narrowPeak2df(narrowPeak_file):
    """
    input a narrowPeak file, output a pandas dataframe with all columns set.
    """
    peak_df = pd.read_table(
        narrowPeak_file,
        names=[
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "signalValue",
            "pValue",
            "qValue",
            "peak",
        ],
        dtype={
            "chrom": str,
            "chromStart": int,
            "chromEnd": int,
            "name": str,
            "score": int,
            "strand": str,
            "signalValue": float,
            "pValue": float,
            "qValue": float,
            "peak": int,
        },
    )
    return peak_df

def broadPeak2df(narrowPeak_file):
    """
    input a narrowPeak file, output a pandas dataframe with all columns set.
    """
    peak_df = pd.read_table(
        narrowPeak_file,
        names=[
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "signalValue",
            "pValue",
            "qValue",
        ],
        dtype={
            "chrom": str,
            "chromStart": int,
            "chromEnd": int,
            "name": str,
            "score": int,
            "strand": str,
            "signalValue": float,
            "pValue": float,
            "qValue": float,
        },
    )
    return peak_df

def bed2df(bed_file):
    """
    input a bed file, output a pandas dataframe with all columns set.
    """
    peak_df = pd.read_table(
        bed_file,
        names=[
            "chrom",
            "chromStart",
            "chromEnd",
        ],
        dtype={
            "chrom": str,
            "chromStart": int,
            "chromEnd": int,
        },
    )
    return peak_df

def annotated_bed2df(bed_file):
    """
    input a bed file, output a pandas dataframe with all columns set.
    """
    peak_df = pd.read_table(
        bed_file,
        names=[
            "chrom",
            "chromStart",
            "chromEnd",
            "annotation"
        ],
        dtype={
            "chrom": str,
            "chromStart": int,
            "chromEnd": int,
            "annotation":str,
        },
    )
    return peak_df


def make_autopct(values):
    """
    prettier pie chart inc. percentages of each slice.
    """
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

def genome_TSS_annotation_prom(genome_path_gtf, TSS_gtf, promlength):
    """
    Takes a ensamble genome gtf file as input, and a specific size for the promoter length
    and returns a smaller gtf file containing only the genes,
    their TSS locations, their strandinfo, and the coordinates of the promoter region.
    """
    Path(f"{TSS_gtf}").touch()
    tempfile = os.getcwd() + "/temp_TSS_genome_annotation.txt"
    Path(tempfile).touch()
    os.remove(f"{TSS_gtf}")
    # i = 0
    with open(f"{tempfile}", "a") as output_gtf:
        with open(genome_path_gtf) as f:
            for _ in range(5):
                next(f)
            for line in f:
                line = line.strip()
                line = line.split()
                gene_name = line[13]
                coding_stuff = line[17]
                strand_info = line[6]
                #             if i == 1000:
                #                 break
                #             i = i + 1
                gene_name = gene_name.rstrip(";")
                coding_stuff = coding_stuff.rstrip(";")
                # coding_stuff = coding_stuff.rstrip(";")
                if str(coding_stuff) == '"protein_coding"':
                    if line[6] == "+":
                        start = line[3]
                        end = int(start) + 1
                        promstart = int(start) - promlength
                        promend = start
                        output_gtf.write(
                            line[0]
                            + "\t"
                            + str(start)
                            + "\t"
                            + str(end)
                            + "\t"
                            + str(promstart)
                            + "\t"
                            + str(promend)
                            + "\t"
                            + strand_info
                            + "\t"
                            + gene_name
                            + "\n"
                        )
                    else:
                        start = int(line[4]) - 1
                        end = line[4]
                        promend = int(start) + promlength
                        promstart = start
                        output_gtf.write(
                            line[0]
                            + "\t"
                            + str(start)
                            + "\t"
                            + str(end)
                            + "\t"
                            + str(promstart)
                            + "\t"
                            + str(promend)
                            + "\t"
                            + strand_info
                            + "\t"
                            + gene_name
                            + "\n"
                        )
    sp.check_call(
        f"nice -5 bedtools sort " f"-i {tempfile} " f"> {TSS_gtf}", shell=True,
    )
    os.remove(tempfile)

def genome_TSS_annotation(genome_path_gtf, TSS_gtf):
    """
    Takes a ensamble genome gtf file as input, and returns a smaller gtf file containing only the genes and
    their TSS locations
    """

    Path(f"{TSS_gtf}").touch()
    tempfile = os.getcwd() + "/temp_TSS_genome_annotation.txt"
    Path(tempfile).touch()
    os.remove(f"{TSS_gtf}")
    # i = 0
    with open(f"{tempfile}", "a") as output_gtf:
        with open(genome_path_gtf) as f:
            for _ in range(5):
                next(f)
            for line in f:
                line = line.strip()
                line = line.split()
                gene_name = line[13]
                coding_stuff = line[17]
                strand_info = line[6]
                #             if i == 1000:
                #                 break
                #             i = i + 1
                gene_name = gene_name.rstrip(";")
                coding_stuff = coding_stuff.rstrip(";")
                # coding_stuff = coding_stuff.rstrip(";")
                if str(coding_stuff) == '"protein_coding"':
                    if line[6] == "+":
                        start = line[3]
                        end = int(start) + 1
                        output_gtf.write(
                            line[0]
                            + "\t"
                            + str(start)
                            + "\t"
                            + str(end)
                            + "\t"
                            + gene_name
                            + "\n"
                        )
                    else:
                        start = int(line[4]) - 1
                        end = line[4]
                        output_gtf.write(
                            line[0]
                            + "\t"
                            + str(start)
                            + "\t"
                            + str(end)
                            + "\t"
                            + gene_name
                            + "\n"
                        )
    sp.check_call(
        f"nice -5 bedtools sort " f"-i {tempfile} " f"> {TSS_gtf}", shell=True,
    )
    os.remove(tempfile)                        




                        
def bedtool_closest(genome_TSS_gtf, bedregion_file):
    """
    Takes a genome TSS gtf file as input (see the 'genome_TSS_GTF' script) together with a bed file with regions. It wil
    subsequently map the regions from the bed file to the closest TSS region
    """
    sp.check_call(
    f"nice -14 bedtools closest "
    f"-k 1 "  # report the clossest gene
    f"-d "  # reports the distance to the gene
    f"-b {genome_TSS_gtf} "
    f"-a {bedregion_file}  "
    f"> tmp_bedtool_closest.txt",
    shell=True)
        
    pd_table = pd.read_table(
    "tmp_bedtool_closest.txt",
    index_col=None,
    names=["chr", "start", "end", "gene_chr", "gene_start", "gene_end", "gene", "region-gene-dist"],
    dtype={
        "chr_e": str,
        "start_e": int,
        "end_e": int,
        "chr": str,
        "start": int,
        "end": int,
        "gene": str,
        "enh-gene-dist": str,
    })
        
    os.remove("tmp_bedtool_closest.txt")
        
    return(pd_table)

def prepend_line(file_name, line):
    """ Insert given string as a new line at the beginning of a file """
    # define name of temporary dummy file
    dummy_file = file_name + ".bak"
    # open original file in read mode and dummy file in write mode
    with open(file_name, "r") as read_obj, open(dummy_file, "w") as write_obj:
        # Write given line to the dummy file
        write_obj.write(line + "\n")
        # Read lines from original file one by one and append them to the dummy file
        for line in read_obj:
            write_obj.write(line)
    # remove original file
    os.remove(file_name)
    # Rename dummy file as the original file
    os.rename(dummy_file, file_name)
                        
def find_active_and_variable_regions(
    old_df, activity_cutoff, variance_cutoff, annotation_string_1, annotation_string_2
):
    """
    first removes all rows with a rowmax below the activity cutoff. 
    Next calculates row variance and keeps rows above the variance cutoff.
    Finally looks if the first or second column is higher, and generates a annotation column.
    if the first column vallue is higher, annotation string 1 is in the annotation column.
    if the second column vallue, annotation string 2 is used in the annotation column.
    
    df = dataframe of activity measurements. with ensemble location (e.g. chr1:100-2000) as a index
    activity_cutoff is a cutoff vallue to specify active regions
    variance_cutoff is a cutoff vallue to specify the variable regions.
    annotation_string_1: outputted if the first column is higher than the second
    annotation_string_2: outputted if the second column is higher than the first
    """
    df = old_df.copy(deep=True)
    df.loc[:, "active"] = df.max(axis=1) > activity_cutoff
    df.loc[:, "var"] = df.iloc[:, [0, 1]].var(axis=1)
    df.loc[:, "variable"] = df.loc[:, "var"] > variance_cutoff
    df.loc[:, "ensmbl_loc"] = df.index

    conditions = [
        (df.loc[:, "active"] == False),
        ((df.loc[:, "active"] == True) & (df.loc[:, "variable"] == False)),
        (
            (df.loc[:, "active"] == True)
            & (df.loc[:, "variable"] == True)
            & (df.iloc[:, 0] > df.iloc[:, 1])
        ),
        (
            (df.loc[:, "active"] == True)
            & (df.loc[:, "variable"] == True)
            & (df.iloc[:, 0] < df.iloc[:, 1])
        ),
    ]

    choices = ["in-active", "common", annotation_string_1, annotation_string_2]

    df["annotation"] = np.select(conditions, choices, default="in-active")

    return df

def TSS_window_to_region(genome_TSS_gtf, bedregion_file, TSS_window):
    """
    Takes a genome TSS gtf file or bedfile as input (see the 'genome_TSS_GTF' script) together with a bed file with regions and
    a window length. 
    It wil subsequently take a window around the TSS, and list all regions overlapping with this window. It will
    furhtermore calculate the distance of te bedregion to the TSS. 
    It will return a dataframe containing all these results.
    """
    sp.check_call(
        f"nice -14 bedtools window "
        f"-w {TSS_window} "
        f"-a {genome_TSS_gtf} "
        f"-b {bedregion_file}  "
        f"> tmp_bedtool_window.txt",
        shell=True,
    )

    pd_table = pd.read_table(
        "tmp_bedtool_window.txt",
        index_col=None,
        names=[
            "chr_TSS",
            "start_TSS",
            "end_TSS",
            "gene_name",
            "Chrom",
            "ChromStart",
            "ChromEnd",
            #"region_annotation",
        ],
        dtype={
            "chr_TSS": str,
            "start_TSS": int,
            "end_TSS": int,
            "gene_name": str,
            "Chrom": str,
            "ChromStart": int,
            "ChromEnd": int,
            "region_annotation": str,
        },
    )

    os.remove("tmp_bedtool_window.txt")

    distance = (
        []
    )  # calculate the distance based on the TSS start end region start site, check which one is upstream for distance calc.
    for x, y in zip(pd_table["start_TSS"], pd_table["ChromStart"]):
        if x > y:
            distance.append(x - y)
        elif y > x:
            distance.append(y - x)
        else:
            distance.append("NaN")

    pd_table["distance"] = distance
    return pd_table


def bedtool_closest_annotated(genome_TSS_gtf, bedregion_file):
    """
    Takes a genome TSS gtf file as input (see the 'genome_TSS_GTF' script) together with a bed file with a 4th annotation column and
    with regions. It wil
    subsequently map the regions from the bed file to the closest TSS region
    """
    sp.check_call(
        f"nice -14 bedtools closest "
        f"-k 1 "  # report the clossest gene
        f"-d "  # reports the distance to the gene
        f"-b {genome_TSS_gtf} "
        f"-a {bedregion_file}  "
        f"> tmp_bedtool_closest.txt",
        shell=True,
    )

    pd_table = pd.read_table(
        "tmp_bedtool_closest.txt",
        index_col=None,
        names=[
            "chr",
            "start",
            "end",
            "annotation",
            "gene_chr",
            "gene_start",
            "gene_end",
            "gene",
            "region-gene-dist",
        ],
        dtype={
            "chr_e": str,
            "start_e": int,
            "end_e": int,
            "chr": str,
            "start": int,
            "end": int,
            "annotation": str,
            "gene": str,
            "enh-gene-dist": str,
        },
    )

    os.remove("tmp_bedtool_closest.txt")

    return pd_table


def distance_weight_region_average(df, weight_dict):
    TSS_dict = {}
    df = df.dropna()
    for gene, KC_int, LSC_int, DIST in zip(
        df["gene_name"], df["KC"], df["LSC"], df["distance"],
    ):
        if gene in TSS_dict:
            old_values = TSS_dict.get(gene)
            new_values = [
                KC_int * weight_dict.get(DIST, 0),
                LSC_int * weight_dict.get(DIST, 0),
                1,
            ]
            # print(old_values)
            # print(new_values)
            summed_values = [
                old_values[i] + new_values[i] for i in range(len(old_values))
            ]
            # print(summed_values)
            TSS_dict.update({gene: summed_values})
        else:
            TSS_dict.update(
                {
                    gene: [
                        KC_int * weight_dict.get(DIST, 0),
                        LSC_int * weight_dict.get(DIST, 0),
                        1,
                    ]
                }
            )

    average_enh_TSS_dict = {}

    for gene in TSS_dict:
        average_enh_TSS_dict.update(
            {
                gene: [
                    (TSS_dict.get(gene)[0] / (TSS_dict.get(gene)[2])),
                    (TSS_dict.get(gene)[1] / (TSS_dict.get(gene)[2])),
                    (TSS_dict.get(gene)[2]),
                ]
            }
        )

    average_enh_df = pd.DataFrame.from_dict(average_enh_TSS_dict, orient="index")
    average_enh_df["mean_int"] = np.mean(average_enh_df.iloc[:, [0, 1]], axis=1)
    #     average_enh_df["FC"] = abs(
    #         np.log2(average_enh_df[0].astype(float))  # / average_enh_df[1])
    #     ).astype(float)
    average_enh_df["FC"] = np.log2(
        average_enh_df[0].astype(float) / average_enh_df[1].astype(float)
    ).astype(float)

    average_enh_df["abs_FC"] = abs(average_enh_df["FC"])

    return average_enh_df


def distance_weight_promoter_region_average(df):
    TSS_dict = {}
    df = df.dropna()

    for gene, KC_int, LSC_int, DIST in zip(
        df["gene_name"], df["KC"], df["LSC"], df["distance"],
    ):
        if gene in TSS_dict:
            old_values = TSS_dict.get(gene)
            new_values = [
                KC_int * 1,
                LSC_int * 1,
                1,
            ]
            # print(old_values)
            # print(new_values)
            summed_values = [
                old_values[i] + new_values[i] for i in range(len(old_values))
            ]
            # print(summed_values)
            TSS_dict.update({gene: summed_values})
        else:
            TSS_dict.update({gene: [KC_int * 1, LSC_int * 1, 1,]})

    average_enh_TSS_dict = {}

    for gene in TSS_dict:
        average_enh_TSS_dict.update(
            {
                gene: [
                    (TSS_dict.get(gene)[0] / (TSS_dict.get(gene)[2])),
                    (TSS_dict.get(gene)[1] / (TSS_dict.get(gene)[2])),
                    (TSS_dict.get(gene)[2]),
                ]
            }
        )

    average_enh_df = pd.DataFrame.from_dict(average_enh_TSS_dict, orient="index")
    average_enh_df["mean_int"] = np.mean(average_enh_df.iloc[:, [0, 1]], axis=1)
    average_enh_df["FC"] = np.log2(
        (average_enh_df[0].astype(float)) / (average_enh_df[1].astype(float))
    ).astype(float)
    average_enh_df["abs_FC"] = np.log2(
        (average_enh_df[0].astype(float)) / (average_enh_df[1].astype(float))
    ).astype(float)

    return average_enh_df

def Ensemble2UCSC(chr_list, db_location, to_ID_type):
    """
    chr_list = "list of chr identifiers, or chr locations (either "1" or "1:1002-3999") that need 
    to be converted (into "chr1" or "chr1:1002-3999")
    db_location = "tab seperated identifier linking file"
    to_ID_type = either 'Ensmble' or UCSC, converts te other to the target ID type.
    
    """
    table = pd.read_table(Ensembl2UCSC_ids)
    if to_ID_type == 'Ensmble':
        table = dict(zip(table.ensembl,table.UCSC))
    elif to_ID_type == 'UCSC':
        table = dict(zip(table.UCSC,table.ensembl))
    else: 
        print('either put Ensmble or UCSC as the to_ID_type input')

    new_ids = []

    if all(chr_list.str.contains(":")):
        print("list contains the ':' symbol, replacing chrom ids in front of :")
        for line in chr_list:
            chrom = line.split(":")[0]
            coordinate = line.split(":")[1]
            chrom = table.get(chrom)
            new_ids.append(chrom + ":" + coordinate)
    else:
        print("no ':' symbol found, asuming only 1 chrom id per line")
        for line in chr_list:
            new_ids.append(table.get(line))

    return new_ids
