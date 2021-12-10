
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



###########################################################################
#### CAST_EiJ -------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=CAST_EiJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/CAST_EiJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64016_190817_211043.subreads.fq.gz $STRAIN_READS/m64016_190820_113556.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

# ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"





###########################################################################
#### CBA_J ----------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=CBA_J
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/CBA_J/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64016_191201_020223.subreads.fq.gz $STRAIN_READS/m64016_200218_143152.subreads.fq.gz $STRAIN_READS/m64016_200221_073449.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"




###########################################################################
#### FVB_NJ ---------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=FVB_NJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/FVB_NJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64089_191108_123542.subreads.fq.gz $STRAIN_READS/m64089_191112_064547.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"



###########################################################################
#### LP_J -----------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=LP_J
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/LP_J/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64089_191113_132427.subreads.fq.gz $STRAIN_READS/m64094_191108_123422.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"




###########################################################################
#### NOD_ShiLtJ -----------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=NOD_ShiLtJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/NOD_ShiLtJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64094_191111_152941.subreads.fq.gz $STRAIN_READS/m64094_191113_060702.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"




###########################################################################
#### NZO_HlLtJ ------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=NZO_HlLtJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/NZO_HlLtJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64089_191109_035327.subreads.fq.gz $STRAIN_READS/m64089_191112_220503.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"




###########################################################################
#### PWK_PhJ --------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=PWK_PhJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/PWK_PhJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64094_191031_141345.subreads.fq.gz $STRAIN_READS/m64094_191102_092811.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"



###########################################################################
#### SPRET_EiJ ------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=SPRET_EiJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/SPRET_EiJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64089_191103_005043.subreads.fq.gz $STRAIN_READS/m64089_191105_065809.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"




###########################################################################
#### WSB_EiJ --------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=WSB_EiJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/WSB_EiJ/pacbio-wgs


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64089_191102_093132.subreads.fq.gz $STRAIN_READS/m64089_191104_154028.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"




###########################################################################
#### C57BL_6NJ ------------------------------------------------------------
###########################################################################

THREADS=16
BTHREADS=$(( $THREADS + 1 )) # minimap2 require and additional thread for I/O
BMEMORY=42000

HOME=/hps/nobackup/production/mousegenomes/emilio/final
STRAIN=C57BL_6NJ
STRAIN_READS=/nfs/production/mousegenomes/rawdata/InbredStrains/C57BL_6NJ/pacbio-wgs/


REFERENCE=/hps/nobackup/production/mousegenomes/emilio/data/Mus_musculus.GRCm39.dna.toplevel.fa.gz
READS="$STRAIN_READS/m64094_191101_181010.subreads.fq.gz $STRAIN_READS/m64089_191031_141215.subreads.fq.gz"

JOB_DIR=$HOME/lsf_out
JOB_ID="$STRAIN"

JOB_ERR=$JOB_DIR/$JOB_ID-align.err
JOB_OUT=$JOB_DIR/$JOB_ID-align.out

OUTPUT_BAM=$HOME/minimap2/$STRAIN.sorted.bam

COMMAND_ALIGNMENT="minimap2 -R '@RG\tID:$STRAIN\tSM:$STRAIN' --MD -Y -t $THREADS -ax map-pb $REFERENCE $READS | samtools view -bS - | samtools sort -T $STRAIN-align -o $OUTPUT_BAM - "

ls $READS && echo $COMMAND_ALIGNMENT

bsub -J "$JOB_ID" -R "select[type==X86_64 && mem > $BMEMORY] rusage[mem=$BMEMORY] span[hosts=1]" -M$BMEMORY -o $JOB_OUT -e $JOB_ERR -n $BTHREADS "$COMMAND_ALIGNMENT"

