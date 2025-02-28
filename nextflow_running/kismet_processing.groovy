//Nextflow pipeline for processing Navin lab methylation libraries//

//TO DO
/*
//amethyst container 
//add bcl-convert
//add bsbolt

//add amethyst initiation (100kb windows, 50kb windows or supplied bed file)
//add scalebio processing flag
//add scalebio indexes for demultiplexing
//copykit redo reference windows to account for bsconversion mappability
*/


// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.flowcellDir = "/volumes/seq/flowcells/MDA/nextseq2000/2024/250127_RM10xMET_RYExome" //Sequencing run flowcell dir
params.src = "/volumes/USR2/Ryan/projects/10x_MET/scmet_nf_processing/src"
params.ref_index="/volumes/USR2/Ryan/projects/10x_MET/ref/hg38_bsbolt"

params.sequencing_cycles="Y151;I10;U16;Y151" // Treat index 2 as UMI just for counting sake
params.cellranger="/volumes/USR2/Ryan/tools/cellranger-atac-2.1.0/"

//library parameters
params.cell_try="5000" //Based on expected cell count from library generation
params.i7_idx="ACTGGTAGAT" //i7 Index (See i7 Indexes in 10xmet_design tab)

//output
params.outname = "250130_10xMET_231_nftest"
params.outdir = "/volumes/USR2/Ryan/projects/10x_MET/experiments/250130_10xmet_231_nf"

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
		Split out Cell ID for N= ${params.cell_try} cells.
		================================================

""".stripIndent()

// BCL TO FASTQ PIPELINE FOR GENERATING SINGLE-CELL FASTQs
process BCL_TO_FASTQ_INIT { 
	//Generate Undetermined Fastq Files from BCL Files.
    //Count GEM indexes and generate a white list for splitting
	//Assumes Y151;I10;U16;Y151 sequencing cycles unless specified as input parameter
	//Assumed ACTGGTAGAT as i7 index unless specified as input parameter
	//TODO This container should be updated to be in the SIF and not local run
	cpus 50
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"

	input:
		path flowcellDir
	output:
		tuple path("initial_gem_idx.txt"), path(flowcellDir)
    script:
		"""
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
	cpus 50
	label 'amethyst'
	containerOptions "--bind ${params.src}:/src/,${params.cellranger}:/cellranger/"
  	publishDir "${params.outdir}/samplesheet", mode: 'copy', overwrite: true, pattern: "samplesheet_gemidx.csv"

	input:
		tuple path(gem_idx), path(flowcellDir)
	output:
		tuple path("samplesheet_gemidx.csv"), path(flowcellDir)
	script:
	"""
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
	cpus 50
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"

	input:
		tuple path(gem_whitelist),path(flowcellDir)
	output:
		path("*_R1_001.fastq.gz")
		path("*_R2_001.fastq.gz")
    script:
		"""
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
		"""
}

// TRIM, ALIGN, and DEDUPLICATE READS
process ADAPTER_TRIM {
	//TRIM READS OF ADAPTERS AND KNOWN METHYLATED REGIONS (GAP FILLS)
	publishDir "${params.outdir}/reports/adapter_trim", mode: 'copy', overwrite: true, pattern: "*.log"
	containerOptions "--bind ${params.src}:/src/,${params.outdir},/volumes/USR2/Ryan/.local/bin:/tools/"
	//TODO This container should be updated to be in the SIF and not local run
	label 'amethyst'

	input:
		tuple path(read1),path(read2)
	output:
		tuple path("*.R1_001.trim.fastq.gz"), path("*.R2_001.trim.fastq.gz")
		//path("*.trim_report.log"), emit: trim_log
	script:
		def sample_name = read1.baseName
		"""
		/tools/cutadapt \\
		-j 1 \
		-a AGATCGGAAGAGCACAC -A CTGTCTCTTATACACAT \
		-U 10 -u 10 \\
		-o ${sample_name}.R1_001.trim.fastq.gz \\
		-p ${sample_name}.R2_001.trim.fastq.gz \\
		$read1 $read2 >> ${sample_name}.trim_report.log 2>> ${sample_name}.trim_report.log
		"""
}

process ALIGN_BSBOLT {
	//ALIGN TRIMMED READS PER CELL
	publishDir "${params.outdir}/reports/alignment", mode: 'copy', overwrite: true, pattern: "*.log"
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"
	label 'amethyst'
	//TODO This container should be updated to be in the SIF and not local run

	input:
		tuple val(sample_name),path(read1),path(read2)
	output:
		tuple val("$sample_name"),path("*.bam")
		//path("*.bsbolt.log"), emit: bsbolt_log
	script:
		"""
		bsbolt Align \\
		-F1 $read1 \\
		-F2 $read2 \\
		-t 1 -OT 1 \\
		-UN -j \\
		-O ${sample_name} \\
		-DB ${params.ref_index} >> ${sample_name}.bsbolt.log 2>> ${sample_name}.bsbolt.log
		"""
}

process MARK_DUPLICATES {
	//MARK DUPLICATE ALIGNMENTS
	//TODO This container should be updated to be in the SIF and not local run
	publishDir "${params.outdir}/reports/markduplicates", mode: 'copy', overwrite: true, pattern: "*.log"
	publishDir "${params.outdir}/sc_bam", mode: 'copy', overwrite: true, pattern: "*.bbrd.bam"
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"
	label 'amethyst'

	input:
		tuple val(sample_name),path(bam)
	output:
		tuple val("$sample_name"),path("*bbrd.bam")
		//path("*markdup.log"), emit: markdup_log
	script:
		"""
		samtools sort -m 10G -n $bam | \\
		samtools fixmate -p -m - - | \\
		samtools sort -m 10G | \\
		samtools markdup --mode s -r -S -s -f ${sample_name}.markdup.log - ${sample_name}.bbrd.bam
	"""
}

process METHYLATION_CALL {
	//CALL CG METHYLATION
	//Split bam file by read names
	publishDir "${params.outdir}/reports/metcalls", mode: 'copy', overwrite: true, pattern: "*.log"
	publishDir "${params.outdir}/sc_metcalls", mode: 'copy', overwrite: true, pattern: "*.h5.gz"

	containerOptions "--bind ${params.src}:/src/,${params.outdir}"
	label 'amethyst'

	input:
		tuple val(sample_name), path(bam)
	output:
		tuple val("$sample_name"),path("*sam")
		//path("*.metcall.log"), emit: metcall_log

	script:
	"""
    samtools index $bam
    bsbolt CallMethylation \\
    -I $bam \\
    -O $sample_name \\
    -ignore-ov -verbose \\
    -min 1 -t 1 -CG \\
    -DB ${params.ref_index} >> ${sample_name}.bsbolt.metcall.log 2>> ${sample_name}.bsbolt.metcall.log

	python /src/premethyst_cgmap_to_h5.py \\
	--input ${sample_name}.CG.map.gz

	"""
}

// CNV PROFILING 
process CNV_CLONES {
	//COPYKIT FOR CLONE CALLING BY CNVS
	//Run CopyKit and output list of bam files by clones
	label 'cnv'
	publishDir "${params.outdir}/cnv_calling", mode: 'copy', pattern: "*{tsv,rds}"
	publishDir "${params.outdir}/plots/cnv", mode: 'copy', pattern: "*pdf"
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"
	cpus 50

	input:
		path bams
	output:
		path("*.scCNA.tsv")
	script:
		"""
		Rscript /src/copykit_cnv_clones.nf.R \\
		--input_dir . \\
		--output_prefix ${params.outname} \\
		--task_cpus ${task.cpus}
		"""
}


process AMETHYST_PROCESSING {
	//INITIATE AMETHYST OBJECT
	//SET H5 LOCATIONS TO OUTPUT DIRECTORY BECAUSE TEMPORARY WORK DIRECTORY IS NOT PERMANENT
	//Split bam file by read names
	publishDir "${params.outdir}/reports/metcalls", mode: 'copy', overwrite: true, pattern: "*.log"
	publishDir "${params.outdir}/sc_metcalls", mode: 'copy', overwrite: true, pattern: "*.h5.gz"

	containerOptions "--bind ${params.src}:/src/,${params.outdir}"
	label 'amethyst'
    
	input:
		tuple val(sample_name), path(bam)
	output:
		tuple val(${sample_name})
		path("*.metcall.log"), emit: metcall_log

	script:
	"""
    samtools index $bam
    bsbolt CallMethylation \\
    -I $bam \\
    -O $sample_name \\
    -ignore-ov -verbose \\
    -min 1 -t 1 -CG \\
    -DB ${params.ref_index} >> ${sample_name}.bsbolt.metcall.log 2>> ${sample_name}.bsbolt.metcall.log

	python /src/premethyst_cgmap_to_h5.py \\
	--input ${sample_name}.CG.map.gz

	"""
}


workflow {
	// BCL TO FASTQ PIPELINE FOR SPLITTING FASTQS
		fq1,fq2 =
		Channel.fromPath(params.flowcellDir) \
		| BCL_TO_FASTQ_INIT \
		| GENERATE_GEM_WHITELIST \
		| BCL_TO_FASTQ_ON_WHITELIST

	//Pair fq1 and fq2 output
        fq1.map { tuple( it.simpleName, it ) }.set{fq1}
		
		fq2.map { tuple( it.simpleName, it ) }.set{fq2}

		fqs = fq1.combine( fq2, by: 0 ).view()

		ADAPTER_TRIM(fqs) 
		| ALIGN_BSBOLT \
		| MARK_DUPLICATES

/*

	//CNV CLONE CALLING
		clone_metadata =  | CNV_CLONES

	//AMETHYST CLONE CALLING
	//METHYLTREE CLONE CALLING
*/
}

/*
example run
source activate #to use more recent version of java

cd /volumes/USR2/Ryan/projects/10x_MET #move to project directory
git clone https://github.com/mulqueenr/scmet_nf_processing #pull github repo

nextflow ./scmet_nf_processing/nextflow_running/kismet_processing.groovy \
-resume \
-with-report \
--flowcellDir /volumes/seq/flowcells/MDA/nextseq2000/2024/250127_RM10xMET_RYExome \
--outname 250130_10xMET_231_nftest

*/


