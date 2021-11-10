import glob
import os
from pathlib import Path

import pandas as pd

bcf_tools = "/home/egarcia/appdir/bcftools/bin/bcftools"
bed_tools = "bedtools"
survivor = "/home/egarcia/workspace/github/SURVIVOR/Debug/SURVIVOR"

base_dir = '/home/egarcia/workspace/github/sv-code/pipeline'

input_calls = glob.glob(base_dir + '/input/calls/*_in.vcf')


# This function always returns an array with: base_dir, extension, part1, part2, part3 ..., part_n
def path_info(file_path, splitter="_"):
    split = os.path.split(file_path)
    splitext = os.path.splitext(split[1])
    name_parts = splitext[0].split(splitter)
    return [split[0], splitext[1]] + name_parts


def prepare_support_data(vcf_file_path, include, target_type, window):
    support_command_sniffles = """
        echo '#CALLER\\tIDX\\tCHROM\\tPOS\\tEND\\tSVLEN\\tSVTYPE\\tFILTER\\tRE\\tAF\\tGT\\tDR\\tDV\\tSVLEN_ABS\\tAF_BIN\\tSVLEN_BIN' > {output} && 
        {bcf_bin} query -i '{include}' -Hf'%CHROM\\t%POS0\\t%END0\\t[%SVLEN]\\t[%SVTYPE]\\t%FILTER\\t[%RE]\\t[%AF]\\t[%GT]\\t[%DR]\\t[%DV]\\n' {vcf_file} | 
        awk -F'\\t' 'BEGIN {{OFS = FS}} $1 ~/^[1-9]*$|^X$/{{ 
            abs=$4<0?-$4:$4; af=(5-(int($7*100)%5)+int($7*100))/100; idx=$1":"$2"-"$3; 
            len_bin=abs<30?"0-30":abs<50?"30-50":abs<100?"50-100":abs<500?"100-500":abs<2000?"500-2000":abs<10000?"2000-10000":abs<50000?"10000-50000":"50000+";  
            print "{caller}",idx,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,abs,af,len_bin
        }}' >> {output}
    """

    support_command_pbsv = """
        echo '#CALLER\\tIDX\\tCHROM\\tPOS\\tEND\\tSVLEN\\tSVTYPE\\tFILTER\\tGT\\tAD_0\\tAD_1\\tDP\\tSVLEN_ABS\\tSVLEN_BIN' > {output} && 
        {bcf_bin} query -i '{include}' -Hf'%CHROM\\t%POS0\\t%END0\\t[%SVLEN]\\t[%SVTYPE]\\t%FILTER\\t[%GT]\\t[%AD{{0}}]\\t[%AD{{1}}]\\t[%DP]\\n' {vcf_file} | 
        awk -F'\\t' 'BEGIN {{OFS = FS}} $1 ~/^[1-9]*$|^X$/{{
            abs=$4<0?-$4:$4; af=(5-(int($7*100)%5)+int($7*100))/100; idx=$1":"$2"-"$3;
            len_bin=abs<30?"0-30":abs<50?"30-50":abs<100?"50-100":abs<500?"100-500":abs<2000?"500-2000":abs<10000?"2000-10000":abs<50000?"10000-50000":"50000+";
            print "{caller}",idx,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,abs,len_bin
        }}' >> {output}
    """

    info = path_info(vcf_file_path)
    strain = info[2].upper()
    caller = info[3].upper()

    output_file_support = base_dir + "/work/merge/raw/" + strain + "_" + caller + "_" + target_type + "_CALLS.tsv"

    if caller == "PBSV":
        command = support_command_pbsv.format(bcf_bin=bcf_tools, include=include, vcf_file=vcf_file_path,
                                              output=output_file_support, caller=caller)
    else:
        command = support_command_sniffles.format(bcf_bin=bcf_tools, include=include, vcf_file=vcf_file_path,
                                                  output=output_file_support, caller=caller)

    print("Generating Support File => " + command)
    stream = os.popen(command)
    print(stream.read())


def prepare_filtered_bed(vcf_file_path, include, target_type, out_dir, suffix="_SRC", window=0):
    bed_command = """
        {bcf_bin} query -i '{include}' -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\n' {vcf_file} | 
         awk -F'\\t' 'BEGIN {{OFS = FS}} $1 ~/^[1-9]*$|^X$/{{print $1,$2,$3-{window},$4+{window}}}' > {output}
    """
    info = path_info(vcf_file_path)
    strain = info[2].upper()
    caller = info[3].upper()

    output_file_bed = base_dir + "/" + out_dir + strain + "_" + caller + "_" + target_type + suffix + ".bed"
    command = bed_command.format(bcf_bin=bcf_tools, include=include, vcf_file=vcf_file_path, window=window,
                                 output=output_file_bed)
    print("Generating BED File => " + command)
    stream = os.popen(command)
    print(stream.read())


def generate_simple_vcf_from_bed():
    bed_to_vcf = "{survivor_bin} bedtovcf {input_bed} {sv_type} {output_file}"
    callers = []
    for input_bed in glob.glob(base_dir + '/work/merge/bed/*_SRC.bed'):
        info = path_info(input_bed)
        b_strain = info[2].upper()
        b_caller = info[3].upper()
        b_type = info[4].upper()

        callers.append(b_strain + "_" + b_caller)
        out_file = base_dir + '/work/merge/vcf/' + b_strain + "_" + b_caller + "_" + b_type + ".vcf"

        command = bed_to_vcf.format(survivor_bin=survivor, input_bed=input_bed, sv_type=b_type, output_file=out_file)

        print("Generating VCF with Survivor => " + command)
        stream = os.popen(command)
        print(stream.read())

    callers = list(dict.fromkeys(callers))
    print(callers)
    view_with_header = """
        {bcf_tools} view {file_in} > {file_out}
    """

    view_no_header = """
        {bcf_tools} view --no-header {file_in} >> {file_out}
    """

    for caller in callers:
        first = True
        file_out = base_dir + '/work/merge/vcf/JOINT_' + caller + ".vcf"
        for file in glob.glob(base_dir + '/work/merge/vcf/' + caller + '_*.vcf'):
            if first:
                first = False
                command = view_with_header.format(bcf_tools=bcf_tools, file_in=file, file_out=file_out)
            else:
                command = view_no_header.format(bcf_tools=bcf_tools, file_in=file, file_out=file_out)

            print("Adding VCF data to joint file => " + command)
            stream = os.popen(command)
            print(stream.read())


def merge_with_survivor():
    file_list = """
        ls {vcf_dir}/JOINT_*.vcf > {vcf_dir}/survivor_input.txt
    """

    command = file_list.format(vcf_dir=base_dir + "/work/merge/vcf")

    print("Generating VCF LIST for Merging => " + command)
    stream = os.popen(command)
    print(stream.read())

    survivor_merge = """
        {survivor_bin} merge {file_list} 150 1 1 1 0 15 {output_vcf}
    """

    command = survivor_merge.format(survivor_bin=survivor, file_list=base_dir + "/work/merge/vcf/survivor_input.txt",
                                    output_vcf=base_dir + "/results/vcf/survivor_merged_calls.vcf")

    print("Generating VCF LIST for Merging => " + command)
    stream = os.popen(command)
    print(stream.read())


def generate_raw_analysis_data():
    survivor_calls = """
        echo "IDX\\tCHROM\\tPOS\\tEND\\tSVLEN\\tSVTYPE\\tSUPP\\tIDX_SRC1\\tIDX_SRC2" > {output_file} &&
        {bcf_tools} query -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\t%SVTYPE\\t%SUPP[\\t%CO]\n' {input_file} | 
        awk  -F'\\t' 'BEGIN {{OFS = FS}} {{idx=$1":"$2"-"$3; print idx,$1,$2,$3,$4,$5,$6,$7,$8}}' >> {output_file} 
    """

    zero_file = base_dir + '/results/analysis/raw/all_calls.tsv'
\
    command = survivor_calls.format(bcf_tools=bcf_tools, input_file=base_dir + '/results/vcf/survivor_merged_calls.vcf',
                                    output_file=zero_file)

    print("Generating Support File for Analysis after Merging => " + command)
    stream = os.popen(command)
    print(stream.read())

    zero_df = pd.read_csv(zero_file, sep='\t', low_memory=False)
    zero_df["IDX_SRC1"] = zero_df["IDX_SRC1"].apply(
        lambda x: x.replace("_", ":").replace("-" + x.split("_")[0] + ":", "-").split(",")[0])
    zero_df["IDX_SRC2"] = zero_df["IDX_SRC2"].apply(
        lambda x: x.replace("_", ":").replace("-" + x.split("_")[0] + ":", "-").split(",")[0])

    zero_df.to_csv(zero_file, sep='\t', index=False)


# Initialize directories
def initialize(dir_path):
    paths = ['/work/merge/raw', '/work/merge/bed', '/work/merge/vcf',
             '/results/vcf', '/results/bed', '/results/analysis/raw',
             '/results/analysis/bed', '/results/analysis/final']
    for path_str in paths:
        Path(dir_path + path_str).mkdir(parents=True, exist_ok=True)


initialize(base_dir)

for vcf_file in input_calls:
    prepare_filtered_bed(vcf_file, 'SVTYPE="DEL" || SVTYPE="DEL/INV"', "DEL", 'work/merge/bed/')
    prepare_support_data(vcf_file, 'SVTYPE="DEL" || SVTYPE="DEL/INV"', "DEL", 0)
    prepare_filtered_bed(vcf_file, 'SVTYPE="INS"', "INS", 'work/merge/bed/')
    prepare_support_data(vcf_file, 'SVTYPE="INS"', "INS", 0)
    prepare_filtered_bed(vcf_file, 'SVTYPE="INV" || SVTYPE="INVDUP"', "INV", 'work/merge/bed/')
    prepare_support_data(vcf_file, 'SVTYPE="INV" || SVTYPE="INVDUP"', "INV", 0)
    prepare_filtered_bed(vcf_file, 'SVTYPE="DUP" || SVTYPE="DUP/INS" || SVTYPE="CNV"', 'work/merge/bed/', "DUP")
    prepare_support_data(vcf_file, 'SVTYPE="DUP" || SVTYPE="DUP/INS" || SVTYPE="CNV"', "DUP", 0)

generate_simple_vcf_from_bed()
merge_with_survivor()

expected_final_vcf = base_dir + "/results/vcf/survivor_merged_calls.vcf"

window = 20

prepare_filtered_bed(expected_final_vcf, 'SVTYPE="DEL"', "DEL", 'results/bed/')
prepare_filtered_bed(expected_final_vcf, 'SVTYPE="DEL"', "DEL", 'results/analysis/bed/', "_" + str(window), window)
prepare_filtered_bed(expected_final_vcf, 'SVTYPE="INS"', "INS", 'results/bed/')
prepare_filtered_bed(expected_final_vcf, 'SVTYPE="INS"', "INS", 'results/analysis/bed/', "_" + str(window), window)
prepare_filtered_bed(expected_final_vcf, 'SVTYPE="INV"', "INV", 'results/bed/')
prepare_filtered_bed(expected_final_vcf, 'SVTYPE="INV"', "INV", 'results/analysis/bed/', "_" + str(window), window)
prepare_filtered_bed(expected_final_vcf, 'SVTYPE="DUP"', "DUP", 'results/bed/')
prepare_filtered_bed(expected_final_vcf, 'SVTYPE="DUP"', "DUP", 'results/analysis/bed/', "_" + str(window), window)

generate_raw_analysis_data()
