# Initial variables
THREADS=12
BMEMORY=120000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog
JOB_DIR=$HOME/pbsv_out

STRAIN=DBA_2J
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$JOB_DIR/$STRAIN-disc.err
JOB_OUT=$JOB_DIR/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"
# bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"




# Initial variables
THREADS=12
BMEMORY=120000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog
JOB_DIR=$HOME/pbsv_out

STRAIN=JF1_MsJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$JOB_DIR/$STRAIN-disc.err
JOB_OUT=$JOB_DIR/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"
# bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"



# Initial variables
THREADS=12
BMEMORY=120000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog
JOB_DIR=$HOME/pbsv_out

STRAIN=129S1_SvImJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$JOB_DIR/$STRAIN-disc.err
JOB_OUT=$JOB_DIR/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"
# bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"



# Initial variables
THREADS=12
BMEMORY=120000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog

STRAIN=AKR_J
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"
# bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"




# Initial variables
THREADS=12
BMEMORY=120000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog

STRAIN=C3H_HeJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"
# bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"



# Initial variables
THREADS=12
BMEMORY=180000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog

STRAIN=CAST_EiJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"	

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"
# bsub -J "$STRAIN" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $THREADS "$COMMAND_DISCOVER"



# Initial variables
THREADS=12
BMEMORY=180000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog

STRAIN=NOD_ShiLtJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"	

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"




# Initial variables
THREADS=12
BMEMORY=120000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog

STRAIN=BALB_cJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"	

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"



# Initial variables
THREADS=12
BMEMORY=120000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog

STRAIN=WSB_EiJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"	

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"




# Initial variables
THREADS=12
BMEMORY=180000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog

STRAIN=A_J
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"	

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"




# Initial variables
THREADS=12
BMEMORY=180000

HOME=/nfs/production/mousegenomes/user/emilio/sv_catalog

STRAIN=NZO_HlLtJ
ALIGNED_READS=$HOME/minimap2/$STRAIN.sorted.bam

JOB_ERR=$HOME/pbsv_out/$STRAIN-disc.err
JOB_OUT=$HOME/pbsv_out/$STRAIN-disc.out

OUTPUT_SVSIG=$HOME/pbsv/svsig/$STRAIN.svsig.gz

COMMAND_DISCOVER="pbsv discover $ALIGNED_READS $OUTPUT_SVSIG"	

# Discoversubmission
bsub -q "dungeon-rh74" -J "$STRAIN-DG" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -n $THREADS "$COMMAND_DISCOVER"




clear && ls -ltrh ../../minimap2/ && ls -ltrh && bjobs