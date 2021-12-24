THREADS=16
BMEMORY=80000

STRAIN=$1

HOME=/hps/nobackup/production/mousegenomes/emilio/final
INPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

# SV Calling variables
OUTPUT_VCF=$HOME/sniffles/$STRAIN-sniffles.vcf
COMMAND_CALL="sniffles --genotype -m $INPUT_BAM -v $OUTPUT_VCF"

JOB_ERR=$HOME/sniffles_out/$STRAIN-sf.err
JOB_OUT=$HOME/sniffles_out/$STRAIN-sf.out

# Variant Call submission
bsub -J "$STRAIN-sf" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_CALL"
