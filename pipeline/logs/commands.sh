
###########################################################################
#### DBA_2J --------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=DBA_2J
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/DBA_2J/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/r64089e_20210528_093647/hifi_reads.fastq.gz $STRAIN_READS/r64089e_20210604_112138_2_B01/hifi_reads.fastq.gz $STRAIN_READS/r64089e_20210604_112138_3_C01/hifi_reads.fastq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-hifi $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"




###########################################################################
#### JF1_MsJ --------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=JF1_MsJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/JF1_MsJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/r64089e_20210422_174152_2_B01/hifi_reads.fastq.gz $STRAIN_READS/r64089e_20210429_165558_1_A01/hifi_reads.fastq.gz $STRAIN_READS/r64089e_20210429_165558_2_B01/hifi_reads.fastq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-hifi $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"




###########################################################################
#### 129S1_SvImJ ----------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=129S1_SvImJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/129S1_SvImJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64016_191026_035320.subreads.fq.gz $STRAIN_READS/m64016_191028_110025.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"



###########################################################################
#### A_J ------------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=A_J
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/A_J/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64089_191111_152802.subreads.fq.gz $STRAIN_READS/m64094_191112_144913.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"




###########################################################################
#### AKR_J ----------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=AKR_J
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/AKR_J/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64016_191025_123544.subreads.fq.gz $STRAIN_READS/m64089_191024_114048.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"



###########################################################################
#### BALB_cJ --------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=BALB_cJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/BALB_cJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64016_190817_055117.subreads.fq.gz $STRAIN_READS/m64016_190819_130947.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"




###########################################################################
#### C3H_HeJ --------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=C3H_HeJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/C3H_HeJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64016_191024_114056.subreads.fq.gz $STRAIN_READS/m64089_191025_132336.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"



