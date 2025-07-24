#!/bin/bash



while read sample;do
#	sample_id=$(basename "$sample" .PASS.filtered2.vcf.gz)
	if [[ $sample == S* ]]; then
		./OptiTypePipeline.py -i /BiO/Hyein/MUTATION2/RNA/${sample}_1.fastq.gz /BiO/Hyein/MUTATION2/RNA/${sample}_2.fastq.gz -r -o ./${sample} --prefix ${sample}
	else
		./OptiTypePipeline.py -i /BiO/Hyein/MUTATION1/RNA/${sample}_1.fastq.gz /BiO/Hyein/MUTATION1/RNA/${sample}_2.fastq.gz -r -o ./${sample} --prefix ${sample}
	fi
done < list
