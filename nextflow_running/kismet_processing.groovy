//Nextflow pipeline for processing Navin lab spatial curio+multiome libraries//

// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.flowcellDir = "/home/rmulqueen/projects/kismet/seq/250127_RM10xMET_RYExome" //Sequencing run flowcell dir
params.src = "/home/rmulqueen/projects/kismet/tools/scmet_nf_processing/src"
params.ref_index="/home/rmulqueen/ref/hg38_bsbolt" 

params.sequencing_cycles="Y151;I10;U16;Y151" // Treat index 2 as UMI just for counting sake
params.cellranger="/home/rmulqueen/tools/cellranger-atac-2.1.0"
params.max_cpus=50

//library parameters
params.cell_try="5000" //Based on expected cell count from library generation
params.i7_idx="ACTGGTAGAT" //i7 Index (See i7 Indexes in 10xmet_design tab)

//output
params.outname = "250130_10xMET_231_nftest"
params.outdir = "/home/rmulqueen/projects/kismet/data/250130_10xmet_231_nf"

log.info """

		================================================
		             kismet Pipeline v1.0
		================================================
		Flowcell Dir : ${params.flowcellDir}
		Sequencing Cycles: ${params.sequencing_cycles}
		Output Prefix : ${params.outname}
		NF Working Dir : ${workflow.launchDir}
		Output Directory : ${params.outdir}

		i7 Index Library Split : ${params.i7_idx}
		Cellranger ATAC install : ${params.cellranger}
		Split out Cell ID for N = ${params.cell_try} cells.

		Max cpus : ${params.max_cpus}
		================================================

""".stripIndent()

// BCL TO FASTQ PIPELINE FOR GENERATING SINGLE-CELL FASTQs
process BCL_TO_FASTQ_INIT { 
	//Generate Undetermined Fastq Files from BCL Files.
    //Count GEM indexes and generate a white list for splitting
	//Assumes Y151;I10;U16;Y151 sequencing cycles unless specified as input parameter
	//bcl-convert requires write access to "/var/logs/bcl-convert", so we just bind a dummy one
	cpus "${params.max_cpus}"
	containerOptions "--bind ${params.outdir}/logs:/var/log/bcl-convert"	
	label 'amethyst'

	input:
		path flowcellDir
	output:
		tuple path("initial_gem_idx.txt"), path(flowcellDir)
    script:
		"""
		source /container_src/container_bashrc
		
		#Generate samplesheet
        echo '[Settings],' > SampleSheet.csv
        echo 'CreateFastqForIndexReads,1' >> SampleSheet.csv
        echo 'OverrideCycles,${params.sequencing_cycles}' >> SampleSheet.csv
        echo '[Data],' >> SampleSheet.csv
        echo 'Sample_ID,index' >> SampleSheet.csv
        echo '${params.outname},${params.i7_idx}' >> SampleSheet.csv

        #Run initial bcl convert and count gem indexes to determine whitelist for splitting
        task_cpus=\$(expr ${task.cpus} / 3)

        bcl-convert \\
        --bcl-input-directory ${flowcellDir} \\
        --bcl-num-conversion-threads \$task_cpus \\
        --bcl-num-compression-threads \$task_cpus \\
        --bcl-num-decompression-threads \$task_cpus \\
        --sample-sheet SampleSheet.csv \\
        --no-lane-splitting true \\
        --output-directory . \\
        --force

        #Count GEM ids for barcodes to keep
        zcat ${params.outname}_S1_R1_001.fastq.gz | \\
        awk 'OFS="\\t" {if(\$1 ~ /^@/) {split(\$1,a,":");print a[8]}}' | \\
        sort -T . --parallel=${task.cpus} --buffer-size=2G | \\
        uniq -c | sort -k1,1n | awk 'OFS="\\t" {print \$1,\$2}' > initial_gem_idx.txt
		"""
}

process GENERATE_GEM_WHITELIST {
	//Take GEM count output from initial Bcl splitting, 
	//generate a new sample sheet for per cell splitting with bcl-convert
	//NEED TO FIX WHITELIST LOCATION, EITHER COPY OR REQUIRE A POINTER TO WHERE CELLRANGER IS INSTALLED
	cpus "${params.max_cpus}"
	label 'amethyst'
	containerOptions "--bind ${params.src}:/src/,${params.cellranger}:/cellranger/"
  	publishDir "${params.outdir}/samplesheet", mode: 'copy', overwrite: true, pattern: "samplesheet_gemidx.csv"

	input:
		tuple path(gem_idx), path(flowcellDir)
	output:
		tuple path("samplesheet_gemidx.csv"), path(flowcellDir)
	script:
	"""
	source /container_src/container_bashrc

	seq_cycles=\$(echo '${params.sequencing_cycles}' | sed 's/U/I/' ) #convert U to I for final cell output
    #make gem specific samplesheet
    python /src/splitcells_whitelist_generator.py \\
    --i7_idx ${params.i7_idx} \\
    --gem_idx $gem_idx \\
    --prefix ${params.outname} \\
    --gem_cutoff ${params.cell_try} \\
	--sequencing_cycles "\${seq_cycles}" \\
	--outdir .
	"""
}

process BCL_TO_FASTQ_ON_WHITELIST { 
	//Generate cell level Fastq Files from BCL Files and generated white list
	//TODO This container should be updated to be in the SIF and not local run
	cpus "${params.max_cpus}"
	containerOptions "--bind ${params.src}:/src/,${params.outdir},${params.outdir}/logs:/var/log/bcl-convert"
	label 'amethyst'
	input:
		tuple path(gem_whitelist),path(flowcellDir)
	output:
		path("*.fastq.gz")
    script:
		"""
		source /container_src/container_bashrc

        #Run final bcl convert to split fastq out per cell
        task_cpus=\$(expr ${task.cpus} / 3)
        bcl-convert \\
        --bcl-input-directory ${flowcellDir} \\
        --bcl-num-conversion-threads \$task_cpus \\
        --bcl-num-compression-threads \$task_cpus \\
        --bcl-num-decompression-threads \$task_cpus \\
		--bcl-only-matched-reads true \\
        --sample-sheet ${gem_whitelist} \\
        --no-lane-splitting true \\
        --output-directory . \\
        --force

		#rename files so simpleName works better
		for file in *R1*fastq.gz; do mv \"\$file\" \"\${file/_R1_/.R1_}\"; done
		for file in *R2*fastq.gz; do mv \"\$file\" \"\${file/_R2_/.R2_}\"; done

		"""
}

// TRIM, ALIGN, and DEDUPLICATE READS
process ADAPTER_TRIM {
	//TRIM READS OF ADAPTERS AND KNOWN METHYLATED REGIONS (GAP FILLS)
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/reports/adapter_trim", mode: 'copy', overwrite: true, pattern: "*.log"
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"
	label 'amethyst'

	input:
		tuple val(cellid),path(read1),path(read2)
	output:
		tuple val(cellid),path("*.R1_001.trim.fastq.gz"), path("*.R2_001.trim.fastq.gz"), emit: fqs
		path("*.trim_report.log"), emit: trim_log
	script:
		"""
		source /container_src/container_bashrc

		cutadapt \\
		-j 1 \\
		-a AGATCGGAAGAGCACAC -A CTGTCTCTTATACACAT \\
		-U 10 -u 10 \\
		-o ${cellid}.R1_001.trim.fastq.gz \\
		-p ${cellid}.R2_001.trim.fastq.gz \\
		$read1 \\
		$read2 >> ${cellid}.trim_report.log 2>> ${cellid}.trim_report.log
		"""
}

process ALIGN_BSBOLT {
	//ALIGN TRIMMED READS PER CELL
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/reports/alignment", mode: 'copy', overwrite: true, pattern: "*.log"
	memory '200 GB'
	label 'amethyst'
	containerOptions "--bind ${params.ref_index}:/ref/"

	input:
		tuple val(cellid),path(read1),path(read2)
	output:
		tuple val(cellid),path("*.bam"), emit: bams
		path("*.bsbolt.log"), emit: bsbolt_log
	script:
		"""
		source /container_src/container_bashrc

		PYTHONPATH=/container_src/bsbolt python -m bsbolt Align \\
		-F1 $read1 \\
		-F2 $read2 \\
		-t 1 -OT 1 \\
		-UN -j \\
		-O ${cellid} \\
		-DB /ref/ >> ${cellid}.bsbolt.log 2>> ${cellid}.bsbolt.log
		"""
}

process MARK_DUPLICATES {
	//MARK DUPLICATE ALIGNMENTS
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/reports/markduplicates", mode: 'copy', overwrite: true, pattern: "*.log"
	publishDir "${params.outdir}/sc_bam", mode: 'copy', overwrite: true, pattern: "*.bbrd.bam"
	label 'amethyst'

	input:
		tuple val(cellid),path(bam)
	output:
		tuple val(cellid),path("*bbrd.bam"), emit: dedup_bams
		path("*markdup.log"), emit: markdup_log
	script:
	"""
		source /container_src/container_bashrc

		samtools sort -m 10G -n $bam | \\
		samtools fixmate -p -m - - | \\
		samtools sort -m 10G | \\
		samtools markdup --mode s -r -S -s -f ${cellid}.markdup.log - ${cellid}.bbrd.bam
	"""
}

process METHYLATION_CALL {
	//CALL CG METHYLATION
	//Split bam file by read names
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/sc_metcalls", mode: 'copy', overwrite: true, pattern: "*.h5.gz"
	publishDir "${params.outdir}/reports/metcalls", mode: 'copy', overwrite: true, pattern: "*.log"
	containerOptions "--bind ${params.ref_index}:/ref/,${params.src}:/src/"
	label 'amethyst'

	input:
		tuple val(cellid), path(bam)
	output:
		tuple val(cellid),path("*h5.gz"), emit: cg_met
		path("*.metcall.log"), emit: metcall_log

	script:
	"""
	source /container_src/container_bashrc

	samtools index $bam
	
	PYTHONPATH=/container_src/bsbolt python -m bsbolt CallMethylation \\
    -I $bam \\
    -O $cellid \\
    -ignore-ov -verbose \\
    -min 1 -t 1 -CG \\
    -DB /ref/ >> ${cellid}.bsbolt.metcall.log 2>> ${cellid}.bsbolt.metcall.log

	python /src/premethyst_cgmap_to_h5.py \\
	--input ${cellid}.CG.map.gz

	"""
}

// CNV PROFILING 
process CNV_CLONES {
	//COPYKIT FOR CLONE CALLING BY CNVS
	//Run CopyKit and output list of bam files by clones
	cpus "${params.max_cpus}"
	label 'cnv'
	publishDir "${params.outdir}/cnv_calling", mode: 'copy', pattern: "*{tsv,rds}"
	publishDir "${params.outdir}/plots/cnv", mode: 'copy', pattern: "*pdf"
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"

	input:
		tuple val(cellid),path(bam)
	output:
		path("*.scCNA.tsv"), emit: clone_tsv
		path("*.scCNA.rds"), emit: clone_rds
		path("*.pdf"), emit: cnv_plots

	script:
		"""
		Rscript /src/copykit_cnv_clones.nf.R \\
		--input_dir . \\
		--output_prefix ${params.outname} \\
		--task_cpus ${task.cpus}
		"""
}

/*
process AMETHYST_PROCESSING {
	//INITIATE AMETHYST OBJECT
	//SET H5 LOCATIONS TO OUTPUT DIRECTORY BECAUSE TEMPORARY WORK DIRECTORY IS NOT PERMANENT
	//Split bam file by read names
	//publishDir "${params.outdir}/reports/metcalls", mode: 'copy', overwrite: true, pattern: "*.log"
	publishDir "${params.outdir}/sc_metcalls", mode: 'copy', overwrite: true, pattern: "*.h5.gz"
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"
	label 'amethyst'
    
	input:
		tuple val(cellid), path(bam)
	output:
		tuple val(${cellid})
		path("*.metcall.log"), emit: metcall_log

	script:
	"""
		source /container_src/container_bashrc
	"""
}
*/

workflow {
	// BCL TO FASTQ PIPELINE FOR SPLITTING FASTQS
		sc_bams = 
		Channel.fromPath(params.flowcellDir) \
		| BCL_TO_FASTQ_INIT \
		| GENERATE_GEM_WHITELIST \
		| BCL_TO_FASTQ_ON_WHITELIST \
		| flatten //combine R1 and R2 to output
		| collate(2) \
		| map { a -> tuple(a[0].simpleName, a[0], a[1]) } \
		| ADAPTER_TRIM

	//Alignment
		ADAPTER_TRIM.out.fqs \
		| ALIGN_BSBOLT

	//Mark duplicates
		ALIGN_BSBOLT.out.bams \
		| MARK_DUPLICATES

	//Call CG methylation
		MARK_DUPLICATES.out.dedup_bams \
		| METHYLATION_CALL

	//CNV CLONE CALLING
		MARK_DUPLICATES.out.dedup_bams \
		| collect \
		| CNV_CLONES

/*
	//AMETHYST CLONE CALLING
	//METHYLTREE CLONE CALLING
*/
}

/*
Example run in readme
*/

