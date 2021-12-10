# Initial variables
THREADS=16
BMEMORY=80000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=DBA_2J

REFERENCE=$HOME/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz

ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_DIR=$HOME/pbsv_out

JOB_ERR=$JOB_DIR/$STRAIN-disc.err
JOB_OUT=$JOB_DIR/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"



# Initial variables
THREADS=12
BMEMORY=160000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=JF1_MsJ

REFERENCE=$HOME/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz

ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_DIR=$HOME/pbsv_out

JOB_ERR=$JOB_DIR/$STRAIN-disc.err
JOB_OUT=$JOB_DIR/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"



# Initial variables
THREADS=8
BMEMORY=160000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=129S1_SvImJ

REFERENCE=$HOME/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz

ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_DIR=$HOME/pbsv_out

JOB_ERR=$JOB_DIR/$STRAIN-disc.err
JOB_OUT=$JOB_DIR/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"
