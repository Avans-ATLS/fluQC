class DashboardText:
    introduction: str = """
FluQC is a dashboard designed for comprehensive quality control of Influenza sequencing data generated by long-read sequencers.
Current QC implemenations:
**Full dataset:**    
    Rough estimation of the percentage of differential interfering particles.
    Mapping statistics generated from samtools coverage.
**Per sample:**
    Readlength distribution per segment.
    Segment coverage heatmap.


For problems, feedback and feature requests, please refer to the [github repo](https://github.com/AVANS-ALST/fluqc)
"""

    dip_explanation: str = """
Differential Interfering Particles (DIPs) are ...

The percentage of DIPs per segment is estimated using the following process:  
1. for each read, calculate alignment ratio: segment length / alignment length.  
2. count the number of reads > THRESHOLD.  
3. calculate percentage of reads passing the threshold.  

The threshold controls how short an alignment can be. The default threshold of 2 counts an aligned read as a putative DIP when the  
alignment is below half the length of the segment it's aligned to. Setting the threshold higher decreases the lenght of alignment  
at which a read is considered a DIP.  
"""
    mapping_explanation: str = """
text\
text\
"""
