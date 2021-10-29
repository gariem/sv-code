import glob
import os
from pathlib import Path

import pandas as pd

bcf_tools = "/home/egarcia/appdir/bcftools/bin/bcftools"
bed_tools = "bedtools"

base_dir = '/home/egarcia/workspace/github/sv-code/pipeline'

input_calls = glob.glob(base_dir + '/input/calls/*.vcf')
input_validated = glob.glob(base_dir + '/input/validated/*.bed')
input_previous = glob.glob(base_dir + '/input/previous/*_chr.bed')


# This function expects a file name in the form: strand_caller.vcf (i.e: /some/path/dba2j_pbsv.vcf)
def file_details_2(file_path):
    file_name = file_path.split("/")[-1]
    split_name = file_name.split("_")
    first = split_name[0]
    second = split_name[1].replace(".vcf", "").replace(".bed", "")
    return file_name, first, second


# This function expects a file name in the form: strand_type_subtype.bed (i.e: /some/path/dba2j_ins_h6.bed)
def file_details_3(file_path):
    file_name = file_path.split("/")[-1]
    split_name = file_name.split("_")
    first = split_name[0]
    second = split_name[1]
    third = split_name[2].replace(".bed", "")
    return file_name, first, second, third


def prepare_exploratory_data(vcf_file_path):
    # TODO: Can infer the command using a list of columns to be included in the query
    query_sniffles = """ 
        echo '#CHROM\\tPOS\\tEND\\tSVLEN\\tSVTYPE\\tFILTER\\tRE\\tAF\\tGT\\tDR\\tDV\\tSVLEN_ABS\\tAF_BIN\\tSVLEN_BIN' > {output} && 
        {bcf_bin} query -Hf'%CHROM\\t%POS0\\t%END0\\t[%SVLEN]\\t[%SVTYPE]\\t%FILTER\\t[%RE]\\t[%AF]\\t[%GT]\\t[%DR]\\t[%DV]\\n' {vcf_file} | 
        awk -F'\\t' 'BEGIN {{OFS = FS}} $1 ~/^[1-9]*$|^X$/{{ 
            abs=$4<0?-$4:$4; af=(5-(int($7*100)%5)+int($7*100))/100; 
            len_bin=abs<30?"0-30":abs<50?"30-50":abs<100?"50-100":abs<500?"100-500":"500+";  
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,abs,af,len_bin
        }}' >> {output}
    """

    query_pbsv = """
        echo '#CHROM\\tPOS\\tEND\\tSVLEN\\tSVTYPE\\tFILTER\\tGT\\tAD\\tDP\\tSVLEN_ABS\\tSVLEN_BIN' > {output} && 
        {bcf_bin} query -Hf'%CHROM\\t%POS0\\t%END0\\t[%SVLEN]\\t[%SVTYPE]\\t%FILTER\\t[%GT]\\t[%AD]\\t[%DP]\\n' {vcf_file} | 
        awk -F'\\t' 'BEGIN {{OFS = FS}} $1 ~/^[1-9]*$|^X$/{{
            abs=$4<0?-$4:$4; af=(5-(int($7*100)%5)+int($7*100))/100; 
            len_bin=abs<30?"0-30":abs<50?"30-50":abs<100?"50-100":abs<500?"100-500":"500+";
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,abs,len_bin
        }}' >> {output}
    """

    file_name, strain, caller = file_details_2(vcf_file_path)
    print("Preparing stats file for [strain: " + strain + " - caller: " + caller + "]")

    output_file = base_dir + "/work/explore/raw/" + strain + "_" + caller + ".tsv"
    if caller == "pbsv":
        command = query_pbsv.format(bcf_bin=bcf_tools, vcf_file=vcf_file_path, output=output_file)
    else:
        command = query_sniffles.format(bcf_bin=bcf_tools, vcf_file=vcf_file_path, output=output_file)

    print("Generating Exploratory RAW File => " + command)
    stream = os.popen(command)
    print(stream.read())


# This function looks for a exploratory file and generates a bed4 from it
def bed_and_intersect_from_exploratory_file(vcf_file_path, window):
    bed_pbsv = """
        awk -F'\\t' 'BEGIN {{OFS = FS}} {{ print $1==\"#CHROM\"?$1:\"chr\"$1,$2-{window},$3+{window},$10 }}' {input} > {output}
    """
    bed_sniffles = """
        awk -F'\\t' 'BEGIN {{OFS = FS}} {{ print $1==\"#CHROM\"?$1:\"chr\"$1,$2-{window},$3+{window},$12 }}' {input} > {output}
    """

    file_name, strain, caller = file_details_2(vcf_file_path)
    exploratory_expected_path = base_dir + "/work/explore/raw/" + strain + "_" + caller + ".tsv"
    generated_bed = base_dir + "/work/explore/bed/" + strain + "_" + caller + ".bed"

    if caller == "pbsv":
        command = bed_pbsv.format(input=exploratory_expected_path, output=generated_bed, window=window)
    else:
        command = bed_sniffles.format(input=exploratory_expected_path, output=generated_bed, window=window)

    print("Generating Exploratory BED File => " + command)
    stream = os.popen(command)
    print(stream.read())

    # Intersect generated file with validated data
    intersect_output = base_dir + '/work/explore/raw/' + strain + "_" + caller + "_intersect_validated.tsv"
    intersect_with(bed_in=generated_bed, files_with=input_validated, intersect_out=intersect_output, window=window)

    # Intersect generated file with previous catalog
    intersect_output = base_dir + '/work/explore/raw/' + strain + "_" + caller + "_intersect_previous.tsv"
    intersect_with(bed_in=generated_bed, files_with=input_previous, intersect_out=intersect_output, window=window)


# Intersect new BED file with experimentally validated data
def intersect_with(bed_in, files_with, intersect_out, window=0):
    # Clear output file and add header
    with open(intersect_out, 'w') as out_file:
        out_file.write('#CHROM\tPOS\tEND\tSVLEN\tTYPE\tCHROM\tPOS\tEND\tSVLEN\tTYPE\tOVERLAP\n')

    for file in files_with:
        file_name, strain, sv_type = file_details_2(file)
        intersect_beds(file, bed_in, intersect_out, window, sv_type)


# Intersect two BED files with -wo option
def intersect_beds(bed_a, bed_b, output, window=0, a_type='N/A', b_type='N/A'):
    intersect = """
        {bed_tools_bin} intersect -a {bed_a} -b {bed_b} -wo | 
        awk -F'\\t' 'BEGIN {{OFS = FS}} {{print $1,$2,$3,$4,toupper("{a_type}"),substr($5,4),$6+{window},$7-{window},$8,"{b_type}",$9}}' >> {output}
    """

    command = intersect.format(bed_tools_bin=bed_tools, bed_a=bed_a, bed_b=bed_b, window=window, output=output,
                               a_type=a_type, b_type=b_type)

    print("Intersecting => " + command)
    stream = os.popen(command)
    print(stream.read())


def merge_tsv_files():
    print("Creating merged Excel file")
    writer = pd.ExcelWriter(base_dir + '/work/explore/merged.xlsx')  # Arbitrary output name
    for tsv_file in glob.glob(base_dir + '/work/explore/raw/*.tsv'):
        sheet_name = tsv_file.replace(tsv_file.replace(".tsv", "").split("_")[0] + "_", "").replace(".tsv", "")
        print("Processing sheet: " + sheet_name)
        df = pd.read_csv(tsv_file, sep='\t', low_memory=False)
        df.to_excel(writer, sheet_name=sheet_name)
    writer.save()


def initialize(dir_path):
    paths = ['/input/calls', '/input/validated', '/input/previous',
             '/work/explore/raw', '/work/explore/bed',
             '/work/vcf', '/work/bed', '/results']
    for path_str in paths:
        Path(dir_path + path_str).mkdir(parents=True, exist_ok=True)


initialize(base_dir)

for vcf_file in input_calls:
    prepare_exploratory_data(vcf_file)
    bed_and_intersect_from_exploratory_file(vcf_file_path=vcf_file, window=50)

merge_tsv_files()
