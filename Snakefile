with open(config['SAMPLES']) as fp:
    samples = fp.read().splitlines()

with open(config['WT']) as fp:
    WT = fp.read().splitlines()

with open(config['nWT']) as fp:
    nWT = fp.read().splitlines()


rule all:
         input:
            #Prepare samples
            #================
            expand("galore/{sample}_1_val_1.fq.gz", sample = samples),
            expand("galore/{sample}_2_val_2.fq.gz", sample = samples),
            expand("{sample}.sam", sample = samples),
            expand("{sample}.bam", sample = samples), 
            expand("{sample}.sorted.bam", sample =samples),
            expand("{sample}.sorted.rmDup.bam", sample =samples),
            expand("{sample}.bigwig", sample = samples),
            "macs/Nfi_1_peaks.narrowPeak",
            "macs/Nfi_2_peaks.narrowPeak",
            "macs/Nfi_1_summits.bed",
            "macs/Nfi_2_summits.bed",
            "Motif_Nfi_1/seq.autonorm.tsv", 
            "Motif_Nfi_2/seq.autonorm.tsv",             
            "nfi_peaks.narrowPeak",
            "mergedPeaksannotated.csv",
            "mergedPeaksStats.txt",
            "Nfi_KEGGpathways.pdf"  
rule trim: 
       input: 
           r1 = "{sample}_1.fastq",
           r2 = "{sample}_2.fastq"
       output: 
          "galore/{sample}_1_val_1.fq.gz",
          "galore/{sample}_2_val_2.fq.gz"
       shell: 
           """
           mkdir -p galore 
           mkdir -p fastqc 
           trim_galore --gzip --retain_unpaired --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """ 
rule align:
              input:               
                  "galore/{sample}_1_val_1.fq.gz",
                  "galore/{sample}_2_val_2.fq.gz"
              params:
                   index=config['INDEX'],
                   mem = config['MEMORY'],
                   cores = config['CORES']
              output:
                   "{sample}.sam",
                   "{sample}_hist.txt" 
              shell:
                   """
                   bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 700 -p {params.cores} -x {params.index} -1 {input[0]} -2 {input[1]} -S {output[0]}  &> {output[1]}
                   """
rule samTobam:
             input: 
                 "{sample}.sam",
             output: 
                 "{sample}.bam"
             shell: 
                   """
                   samtools view -bS {input} > {output}
                   """


rule sort:
             input:
                  "{sample}.bam"
             output:
                    "{sample}.sorted.bam"
             shell:
                   """
                   picard SortSam I={input}  O={output} SORT_ORDER=coordinate
                   """


rule remove_duplicates:
       input: 
        "{sample}.sorted.bam"
       output: 
         "{sample}.sorted.rmDup.bam",
         "{sample}.rmDup.txt"
       shell: 
           """
            picard MarkDuplicates I={input} O={output[0]} REMOVE_DUPLICATES=true METRICS_FILE={output[1]} 
           """

rule index: 
      input: 
         "{sample}.sorted.rmDup.bam"
      output: 
         "{sample}.sorted.rmDup.bam.bai"
      shell: 
          """
          samtools index {input} 
          """ 
rule bamCoverage: 
       input: 
        "{sample}.sorted.rmDup.bam",
        "{sample}.sorted.rmDup.bam.bai" 
       output: 
        "{sample}.bigwig" 
       params: 
         genome_size = config['Genome_Size'], 
         binsize = config['BINSIZE'], 
         num_processors = config['Num_Processors'] 
       shell: 
          """ 
          bamCoverage -b {input[0]} -p {params.num_processors}  --normalizeUsing RPGC --effectiveGenomeSize {params.genome_size} --binSize {params.binsize} -o {output} 
          """ 

rule macs_bed: 
      input: 
         expand("{sample}.sorted.rmDup.bam", sample =samples)
      output: 
          "macs/Nfi_1_peaks.narrowPeak",
          "macs/Nfi_2_peaks.narrowPeak",
          "macs/Nfi_1_summits.bed",
          "macs/Nfi_2_summits.bed" 
      shell: 
           """
           macs2 callpeak -t SRR15632726.sorted.rmDup.bam -c SRR15632728.sorted.rmDup.bam -f BAMPE -g mm --outdir macs -n Nfi_1
           macs2 callpeak -t SRR15632727.sorted.rmDup.bam -c SRR15632728.sorted.rmDup.bam -f BAMPE -g mm --outdir macs -n Nfi_2
           """


rule annotateNarrowPeaks: 
      input: 
          "Nfi_1_peaks.narrowPeak",
          "Nfi_2_peaks.narrowPeak"
      params: 
           genome= config['GENOME'], 
           gtf = config['GTF']  
      output: 
         "Nfi_1.annotatednarrowpeaks", 
         "Nfi_1.annotatednarrowpeaks.stats",
         "Nfi_2.annotatednarrowpeaks",   
         "Nfi_2.annotatednarrowpeaks.stats",
      shell: 
          """
          annotatePeaks.pl {input[[0]} {params.genome} -gtf {params.gtf}   -annStats {output[1]}  > {output[0]} 
          annotatePeaks.pl {input[[1]} {params.genome} -gtf {params.gtf}   -annStats {output[3]}  > {output[2]}   
          """ 


rule findMotifs:
      input:
          "macs/Nfi_1_summits.bed", 
          "macs/Nfi_2_summits.bed"
      params:
          genome = config['GENOME'],
          output_dir1 = "Motif_1",
          output_dir2 = "Motif_2"
      output:
          "Motif_Nfi_1/seq.autonorm.tsv",
          "Motif_Nfi_2/seq.autonorm.tsv"
      shell:
         """
          findMotifsGenome.pl {input[0]} {params.genome} {params.output_dir1} -size 200 -mask
          findMotifsGenome.pl {input[1]} {params.genome} {params.output_dir2} -size 200 -mask
         """


rule sharedPeaks: 
     input: 
           "macs/Nfi_1_peaks.narrowPeak",
           "macs/Nfi_2_peaks.narrowPeak"
     output: 
           "nfi_peaks.narrowPeak"
     shell: 
       """
        bedtools intersect -a {input[0]}  -b {input[1]}  -wo > {output} 
       """


rule annotateSharedPeaks: 
    input: 
       "nfi_peaks.narrowPeak"        
    output: 
       "mergedPeaksannotated.csv", 
       "mergedPeaksStats.txt" 
    params: 
      genome= config['GENOME'],
      gtf = config['GTF']
    shell: 
      """
       annotatePeaks.pl {input} {params.genome} -gtf {params.gtf}   -annStats {output[1]}  > {output[0]} 
      """


rule pathway: 
    input: 
       "nfi_peaks.narrowPeak"
    output: 
       "Nfi_KEGGpathways.pdf" 
    shell: 
        "Rscript pathways.R "
