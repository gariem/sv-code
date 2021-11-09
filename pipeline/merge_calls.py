import glob
import os
from pathlib import Path

bcf_tools = "bcftools"
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


def prepare_filtered_bed_and_support(vcf_file_path, include, target_type, window):
    bed_command = """
        {bcf_bin} query -i '{include}' -f'%CHROM\\t%POS0\\t%END0\\t[%SVLEN]\n' {vcf_file} | 
         awk -F'\\t' 'BEGIN {{OFS = FS}} $1 ~/^[1-9]*$|^X$/{{print $1,$2,$3-{window},$4+{window}}}' > {output}
    """

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

    output_file_bed = base_dir + "/work/merge/bed/" + strain + "_" + caller + "_" + target_type + "_SRC.bed"
    command = bed_command.format(bcf_bin=bcf_tools, include=include, vcf_file=vcf_file_path, window=0,
                                 output=output_file_bed)
    print("Generating BED File => " + command)
    stream = os.popen(command)
    print(stream.read())

    output_file_bed = base_dir + "/work/merge/bed/" + strain + "_" + caller + "_" + target_type + "_WIDE.bed"
    command = bed_command.format(bcf_bin=bcf_tools, include=include, vcf_file=vcf_file_path, window=window,
                                 output=output_file_bed)
    print("Generating BED File => " + command)
    stream = os.popen(command)
    print(stream.read())

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


def generate_simple_vcf_from_bed():
    bedtovcf_command = "{survivor_bin} bedtovcf {input_bed} DEL {output_file}"

    for input_bed in glob.glob(base_dir + '/work/merge/bed/*_SRC.bed'):
        info = path_info(input_bed)
        b_strain = info[2].upper()
        b_caller = info[3].upper()
        b_type = info[4].upper()

        out_file = base_dir + 'work/merge/vcf/' + b_strain + "_" + b_caller + "_" + b_type + ".vcf"

        command = bedtovcf_command.format(survivor_bin=survivor, input_bed=input_bed, output_file=out_file)

        print("Generating VCF with Survivor => " + command)
        stream = os.popen(command)
        print(stream.read())

# Initialize directories
def initialize(dir_path):
    paths = ['/work/merge/raw', '/work/merge/bed', '/work/merge/vcf']
    for path_str in paths:
        Path(dir_path + path_str).mkdir(parents=True, exist_ok=True)


initialize(base_dir)

# for vcf_file in input_calls:
#     prepare_filtered_bed_and_support(vcf_file, 'SVTYPE="DEL" || SVTYPE="DEL/INV"', "DEL", 20)
#     prepare_filtered_bed_and_support(vcf_file, 'SVTYPE="INS"', "INS", 20)
#     prepare_filtered_bed_and_support(vcf_file, 'SVTYPE="INV" || SVTYPE="INVDUP"', "INV", 20)
#     prepare_filtered_bed_and_support(vcf_file, 'SVTYPE="DUP" || SVTYPE="DUP/INS" || SVTYPE="CNV"', "DUP", 20)

generate_simple_vcf_from_bed()