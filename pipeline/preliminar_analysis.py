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


# This function always returns an array with: base_dir, extension, part1, part2, part3 ..., part_n
def path_info(file_path, splitter="_"):
    split = os.path.split(file_path)
    splitext = os.path.splitext(split[1])
    name_parts = splitext[0].split(splitter)
    return [split[0], splitext[1]] + name_parts


# Generate TSV files with human readable info from VCF
def prepare_exploratory_data(vcf_file_path):
    query_sniffles = """ 
        echo '#CALLER\\tIDX\\tCHROM\\tPOS\\tEND\\tSVLEN\\tSVTYPE\\tFILTER\\tRE\\tAF\\tGT\\tDR\\tDV\\tSVLEN_ABS\\tAF_BIN\\tSVLEN_BIN' > {output} && 
        {bcf_bin} query -Hf'%CHROM\\t%POS0\\t%END0\\t[%SVLEN]\\t[%SVTYPE]\\t%FILTER\\t[%RE]\\t[%AF]\\t[%GT]\\t[%DR]\\t[%DV]\\n' {vcf_file} | 
        awk -F'\\t' 'BEGIN {{OFS = FS}} $1 ~/^[1-9]*$|^X$/{{ 
            abs=$4<0?-$4:$4; af=(5-(int($7*100)%5)+int($7*100))/100; idx=$1":"$2"-"$3; 
            len_bin=abs<30?"0-30":abs<50?"30-50":abs<100?"50-100":abs<500?"100-500":abs<2000?"500-2000":abs<10000?"2000-10000":abs<50000?"10000-50000":"large";  
            print "{caller}",idx,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,abs,af,len_bin
        }}' >> {output}
    """

    query_pbsv = """
        echo '#CALLER\\tIDX\\tCHROM\\tPOS\\tEND\\tSVLEN\\tSVTYPE\\tFILTER\\tGT\\tAD\\tDP\\tSVLEN_ABS\\tSVLEN_BIN' > {output} && 
        {bcf_bin} query -Hf'%CHROM\\t%POS0\\t%END0\\t[%SVLEN]\\t[%SVTYPE]\\t%FILTER\\t[%GT]\\t[%AD]\\t[%DP]\\n' {vcf_file} | 
        awk -F'\\t' 'BEGIN {{OFS = FS}} $1 ~/^[1-9]*$|^X$/{{
            abs=$4<0?-$4:$4; af=(5-(int($7*100)%5)+int($7*100))/100; idx=$1":"$2"-"$3;
            len_bin=abs<30?"0-30":abs<50?"30-50":abs<100?"50-100":abs<500?"100-500":abs<2000?"500-2000":abs<10000?"2000-10000":abs<50000?"10000-50000":"large";
            print "{caller}",idx,$1,$2,$3,$4,$5,$6,$7,$8,$9,abs,len_bin
        }}' >> {output}
    """

    info = path_info(vcf_file_path)
    strain = info[2]
    caller = info[3]
    print("Preparing stats file for [strain: " + strain + " - caller: " + caller + "]")

    output_file = base_dir + "/work/explore/raw/" + strain + "_" + caller + ".tsv"
    if caller == "pbsv":
        command = query_pbsv.format(bcf_bin=bcf_tools, vcf_file=vcf_file_path, output=output_file, caller=caller)
    else:
        command = query_sniffles.format(bcf_bin=bcf_tools, vcf_file=vcf_file_path, output=output_file, caller=caller)

    print("Generating Exploratory RAW File => " + command)
    stream = os.popen(command)
    print(stream.read())


# This function looks for a exploratory file and generates a bed4 from it
def bed_and_intersect_from_exploratory_file(vcf_file_path, window):
    bed_pbsv = """
        awk -F'\\t' 'BEGIN {{OFS = FS}} {{ print $3==\"CHROM\"?"#"$3:\"chr\"$3,$4-{window},$5+{window},$12 }}' {input} > {output}
    """
    bed_sniffles = """
        awk -F'\\t' 'BEGIN {{OFS = FS}} {{ print $3==\"CHROM\"?"#"$3:\"chr\"$3,$4-{window},$5+{window},$14 }}' {input} > {output}
    """

    info = path_info(vcf_file_path)
    strain = info[2]
    caller = info[3]

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
    intersect_with(bed_in=generated_bed, files_with=input_validated, intersect_out=intersect_output, window=window,
                   caller=caller)

    # Intersect generated file with previous catalog
    intersect_output = base_dir + '/work/explore/raw/' + strain + "_" + caller + "_intersect_previous.tsv"
    intersect_with(bed_in=generated_bed, files_with=input_previous, intersect_out=intersect_output, window=window,
                   caller=caller)


# Intersect new BED file with experimentally validated data
def intersect_with(bed_in, files_with, intersect_out, window=0, caller=''):
    # Clear output file and add header
    with open(intersect_out, 'w') as out_file:
        out_file.write(
            '#CALLER\tIDX\tCHROM_A\tPOS_A\tEND_A\tSVLEN_A\tTYPE_A\tCHROM_B\tPOS_B\tEND_B\tSVLEN_B\tTYPE_B\tOVERLAP\n')

    for file in files_with:
        info = path_info(file)
        sv_type = info[3]
        intersect_beds(file, bed_in, intersect_out, window, sv_type, caller=caller)


# Intersect two BED files with -wo option
def intersect_beds(bed_a, bed_b, output, window=0, a_type='', b_type='', caller=''):
    intersect = """
        {bed_tools_bin} intersect -a {bed_a} -b {bed_b} -wo | 
        awk -F'\\t' 'BEGIN {{OFS = FS}} {{
        idx=substr($5,4)":"$6+{window}"-"$7-{window}
        print "{caller}",idx,$1,$2,$3,$4,toupper("{a_type}"),substr($5,4),$6+{window},$7-{window},$8,"{b_type}",$9}}' >> {output}
    """

    command = intersect.format(bed_tools_bin=bed_tools, bed_a=bed_a, bed_b=bed_b, window=window, output=output,
                               a_type=a_type, b_type=b_type, caller=caller)

    print("Intersecting => " + command)
    stream = os.popen(command)
    print(stream.read())


# Combine TSV files into a single XLSX file
def merge_tsv_files():
    print("Creating merged Excel file")
    writer = pd.ExcelWriter(base_dir + '/work/explore/merged.xlsx')  # Arbitrary output name
    for tsv_file in glob.glob(base_dir + '/work/explore/raw/*.tsv'):
        sheet_name = "_".join(path_info(tsv_file)[2:])
        print("Processing sheet: " + sheet_name)
        df = pd.read_csv(tsv_file, sep='\t', low_memory=False)
        df.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()


# Join validated intersections with original new calls
def join_with_intersect(set_identifier):
    for tsv_file in glob.glob(base_dir + '/work/explore/raw/*_intersect_' + set_identifier + '.tsv'):
        info = path_info(tsv_file)
        folder = info[0]
        strain = info[2]
        caller = info[3]

        original_calls = folder + "/" + strain + "_" + caller + ".tsv"
        join_output = folder + "/" + strain + "_" + caller + "_join_" + set_identifier + ".tsv"

        print("Joining new reads with " + set_identifier + " intersect: " + tsv_file + " <--> " + original_calls)

        intersect_df = pd.read_csv(tsv_file, sep='\t', low_memory=False)
        calls_df = pd.read_csv(original_calls, sep='\t', low_memory=False)

        join_df = pd.merge(intersect_df, calls_df, on="IDX")
        join_df.to_csv(join_output, sep='\t', index=False)


# Initialize directories
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

# merge_tsv_files()
join_with_intersect("validated")
join_with_intersect("previous")
