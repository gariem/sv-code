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
THREADS=8
BMEMORY=120000

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
BMEMORY=120000

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



# Initial variables
THREADS=8
BMEMORY=120000
HOME=/hps/nobackup/production/mousegenomes/emilio/final
REFERENCE=$HOME/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz

STRAIN=A_J
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"



# Initial variables
THREADS=8
BMEMORY=180000
HOME=/hps/nobackup/production/mousegenomes/emilio/final
REFERENCE=$HOME/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz

STRAIN=AKR_J
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"



# Initial variables
THREADS=8
BMEMORY=120000
HOME=/hps/nobackup/production/mousegenomes/emilio/final
REFERENCE=$HOME/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz

STRAIN=BALB_cJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"



# Initial variables
THREADS=8
BMEMORY=180000
HOME=/hps/nobackup/production/mousegenomes/emilio/final
REFERENCE=$HOME/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz

STRAIN=C3H_HeJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"




# Initial variables
THREADS=8
BMEMORY=180000
HOME=/hps/nobackup/production/mousegenomes/emilio/final
REFERENCE=$HOME/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz

STRAIN=CAST_EiJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"
