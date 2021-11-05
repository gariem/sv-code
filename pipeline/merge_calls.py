import glob
import os

bcf_tools = "/home/egarcia/appdir/bcftools/bin/bcftools"
bed_tools = "bedtools"

base_dir = '/home/egarcia/workspace/github/sv-code/pipeline'

input_calls = glob.glob(base_dir + '/input/calls/*.vcf')


# This function always returns an array with: base_dir, extension, part1, part2, part3 ..., part_n
def path_info(file_path, splitter="_"):
    split = os.path.split(file_path)
    splitext = os.path.splitext(split[1])
    name_parts = splitext[0].split(splitter)
    return [split[0], splitext[1]] + name_parts


def prepare_filtered_bed(vcf_file_path, include, target_type):
    query_command = """
        {bcf_bin} query -i '{include}' -f'%CHROM\\t%POS0\\t%END0\\t[%SVLEN]\n' {vcf_file} > {output}
    """

    info = path_info(vcf_file_path)
    strain = info[2].upper()
    caller = info[3].upper()

    output_file = base_dir + "/work/bed/" + strain + "_" + caller + "_" + target_type + ".bed"

    command = query_command.format(bcf_bin=bcf_tools, include=include, vcf_file=vcf_file_path, output=output_file)

    print("Generating BED File => " + command)
    stream = os.popen(command)
    print(stream.read())


for vcf_file in input_calls:
    prepare_filtered_bed(vcf_file, 'SVTYPE="INS"', "INS")
    prepare_filtered_bed(vcf_file, 'SVTYPE="INS"', "DEL")
