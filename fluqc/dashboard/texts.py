class DashboardText:
    """Markdown text to use in dashboard"""

    introduction: str = """
FluQC is a dashboard designed for comprehensive quality control of Influenza sequencing data generated by long-read sequencers.  

### Current implementations:
[Summary](/page-summary):\n    
A summary table showing general statistics and subtypes of each sample. \n 
[Mapping Statistics](/page-mapping):\n
Heatmap of mapping statistics per sample, per segment.  \n
[In-depth Sample view](/page-sample):\n  
Bivariate length v.s. quality figure, readlengths per segment and coverage histograms per sample. \n  
[DIPs](/page-dip):\n  
A heatmap showing a (rather rudimentary) estimation of the percentage of Defective 
Interfering Particles (DIPs) in the sample, per segment. \n
   
### Have some feedback?
For problems, feedback and feature requests, please refer to the [github repo](https://github.com/AVANS-ALST/fluqc).  
Feel free to post a new issue!
"""
    epilogue: str = """
FluQC was developed within the [MOEDIG consortium](https://www.sia-projecten.nl/project/moleculaire-epidemiologische-diagnostiek-voor-interventiemogelijkheden-ten-aanzien-van-griep-influenza-op-varkensbedri) 
at research group [Analysis Techniques in the Life Sciences](https://www.avans.nl/onderzoek/expertisecentra/perspectief-in-gezondheid/lectoraten/analysetechnieken-in-de-life-sciences), 
part of [Centre of Expertise Perspective in Health](https://www.avans.nl/onderzoek/expertisecentra/perspectief-in-gezondheid).
"""
    table_text: str = """
## Summary Table
The table below summarizes general statistics of the fastq data given to the dashboard.
The minimum depth is the lowest average depth over a 15bp window across the genome.
"""

    dip_explanation: str = """
Defective Interfering Particles (DIPs) are a phenomenon in viral populations where a virus is missing a part of its genome, 
frequently observed in cultured Influenza viruses. DIPs are thought to be generated by the viral polymerase, which can make errors
during replication. These errors can lead to the generation of truncated viral genomes.

The percentage of DIPs per segment is estimated using the following calculation:  
1. for each read, calculate alignment ratio: segment length / alignment length.  
2. count the number of reads > THRESHOLD.  
3. calculate percentage of reads passing the threshold.  

The threshold controls how short an alignment can be. The default threshold of 2 counts an aligned read as a putative DIP when the  
alignment length is less than half the length of the segment it's aligned to. Setting the threshold higher decreases the length of alignment  
at which a read is considered a DIP.  
"""
    mapping_explanation: str = """   
### Description
This page shows a heatmap of mapping statistics derived from [samtools coverage](https://www.htslib.org/doc/samtools-coverage.html).  
Statistics are shown per sample, per segment, with samples on the x-axis and segments on the y-axis. If a sample has no reads for 
a particular segment, the corresponding cell in the heatmap will be empty.  
  
### Statistics
Choose from one of the following statistics to see the results in the heatmap:  
**numreads**: The number of reads mapped against each segment.  
**covbases**: The number of bases contained within the mapped reads per segment.  
**coverage**: The percentage of each segment covered by reads. (percentage of reference segment length)  
**meandepth**: The average number of reads covering each base in the segment.  
**meanbaseq**: The average quality of reads covering the segment.  
**meanmapq**: The average mapping quality of all mapped reads.  
"""
