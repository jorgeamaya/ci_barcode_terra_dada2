version 1.0

workflow amplicon_decontamination_detect {
	input {
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"
		Int read_maxlength = 200
		Int pairread_minlength = 100
		Int merge_minlength = 100
		Int joined_threshold = 1000
		Float contamination_threshold = 0.5
		File barcodes_file
		File pr1 
		File pr2 
		String Class = "parasite"
		String maxEE = "5,5"
		String trimRight = "0,0"
		Int minLen = 30
		String truncQ = "5,5"
		String matchIDs = "0"
		Int max_consist = 10
		Float omegaA = 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
		String saveRdata = ""
		Int justConcatenate = 0
		Int maxMismatch = 0
		File path_to_DADA2
		File overlap_pr1 
		File overlap_pr2 
		File reference
		String adjust_mode = "absolute"
		File path_to_snv 
		String no_ref = "False"
		File reference2 
		String strain = "3D7"
		String strain2 = "DD2"
		String polyN = "5"
		String min_reads = "0"
		String min_samples = "0"
		String max_snv_dist = "-1"
		String max_indel_dist = "-1"
		String include_failed = "False"
		String exclude_bimeras = "False"
		String amp_mask = "None"
		String verbose = "False"
	}
	call ampseq_bbmerge_process {
		input:
			path_to_fq = path_to_fq,
			path_to_flist = path_to_flist,
			pattern_fw = pattern_fw,
			pattern_rv = pattern_rv,
			read_maxlength = read_maxlength,
			pairread_minlength = pairread_minlength,
			merge_minlength = merge_minlength,
			joined_threshold = joined_threshold,
			contamination_threshold = contamination_threshold,
			barcodes_file = barcodes_file,
			pr1 = pr1,
			pr2 = pr2,
			Class = Class,
			maxEE = maxEE,
			trimRight = trimRight,
			minLen = minLen,
			truncQ = truncQ,
			matchIDs = matchIDs,
			max_consist = max_consist,
			omegaA = omegaA,
			saveRdata = saveRdata,
			justConcatenate = justConcatenate,
			maxMismatch = maxMismatch,
			path_to_DADA2 = path_to_DADA2,
			overlap_pr1 = overlap_pr1,
			overlap_pr2 = overlap_pr2,
			reference = reference,
			adjust_mode = adjust_mode,
			path_to_snv = path_to_snv,
			no_ref = no_ref,
			reference2 = reference2,
			strain = strain,
			strain2 = strain2,
			polyN = polyN,
			min_reads = min_reads,
			min_samples = min_samples,
			max_snv_dist = max_snv_dist,
			max_indel_dist = max_indel_dist,
			include_failed = include_failed,
			exclude_bimeras = exclude_bimeras,
			amp_mask = amp_mask,
			verbose = verbose
	}

	output {
		File rawfilelist_f = ampseq_bbmerge_process.rawfilelist
		File missing_files_f = ampseq_bbmerge_process.missing_files
		File merge_tar_f = ampseq_bbmerge_process.merge_tar
		File bbmergefields_f = ampseq_bbmerge_process.bbmergefields
		File BBmerge_performance_absolute_report_f = ampseq_bbmerge_process.BBmerge_performance_absolute_report
		File BBmerge_performance_percentage_report_f = ampseq_bbmerge_process.BBmerge_performance_percentage_report
		File BBmerge_performace_absolute_discarded_f = ampseq_bbmerge_process.BBmerge_performace_absolute_discarded
		File Barcode_report_abs_f = ampseq_bbmerge_process.Barcode_report_abs
		File Barcode_report_per_f = ampseq_bbmerge_process.Barcode_report_per
		File Insert_size_f = ampseq_bbmerge_process.Insert_size
		File Match_report_abs_f = ampseq_bbmerge_process.Match_report_abs
		File Match_report_per_f = ampseq_bbmerge_process.Match_report_per
		File barcodes_report_bbmerge_f = ampseq_bbmerge_process.barcodes_report_bbmerge
		File hamming_distances_forward_f = ampseq_bbmerge_process.hamming_distances_forward
		File hamming_distances_reverse_f = ampseq_bbmerge_process.hamming_distances_reverse
		#?REMOVE
		#File stdout_string = ampseq_bbmerge_process.ampseq_bbmerge_process_stdout
	}
}

task ampseq_bbmerge_process {
	input {
		String path_to_fq 
		File path_to_flist
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"
		Int read_maxlength = 200
		Int pairread_minlength = 100
		Int merge_minlength = 100
		Int joined_threshold = 1000
		Float contamination_threshold = 0.5
		File barcodes_file
		File pr1 
		File pr2
		String Class = "parasite"
		String maxEE = "5,5"
		String trimRight = "0,0"
		Int minLen = 30
		String truncQ = "5,5"
		String matchIDs = "0"
		Int max_consist = 10
		Float omegaA = 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
		String saveRdata = ""
		Int justConcatenate = 0
		Int maxMismatch = 0
		File path_to_DADA2
		File overlap_pr1 
		File overlap_pr2
		File reference
		String adjust_mode = "absolute"
		File path_to_snv 
		String no_ref = "False"
		File reference2
		String strain = "3D7"
		String strain2 = "DD2"
		String polyN = "5"
		String min_reads = "0"
		String min_samples = "0"
		String max_snv_dist = "-1"
		String max_indel_dist = "-1"
		String include_failed = "False"
		String exclude_bimeras = "False"
		String amp_mask = "None"
		String verbose = "False"
	}

	Map[String, String] in_map = {
		"path_to_fq": "fq_dir",
		"path_to_flist": sub(path_to_flist, "gs://", "/cromwell_root/"),
		"pattern_fw": pattern_fw,
		"pattern_rv": pattern_rv,
		"read_maxlength": read_maxlength,
		"pairread_minlength": pairread_minlength,
		"merge_minlength": merge_minlength,
		"joined_threshold": joined_threshold,
		"contamination_threshold": contamination_threshold,
		"barcodes_file": barcodes_file,
		"pr1": pr1,
		"pr2": pr2,
		"Class": Class,
		"maxEE": maxEE,
		"trimRight": trimRight,
		"minLen": minLen,
		"truncQ": truncQ,
		"matchIDs": matchIDs,
		"max_consist": max_consist,
		"omegaA": omegaA,
		"saveRdata": saveRdata,
		"justConcatenate": justConcatenate,
		"maxMismatch": maxMismatch,
		"path_to_DADA2": path_to_DADA2,
		"overlap_pr1": overlap_pr1,
		"overlap_pr2": overlap_pr2,
		"reference": reference,
		"adjust_mode": adjust_mode,
		"path_to_snv": path_to_snv,
		"no_ref": no_ref,
		"reference2": reference2,
		"strain": strain,
		"strain2": strain2,
		"polyN": polyN,
		"min_reads": min_reads,
		"min_samples": min_samples,
		"max_snv_dist": max_snv_dist,
		"max_indel_dist": max_indel_dist,
		"include_failed": include_failed,
		"exclude_bimeras": exclude_bimeras,
		"amp_mask": amp_mask,
		"verbose": verbose
	}
	File config_json = write_json(in_map)
	command <<<
	set -euxo pipefail
	#set -x
	mkdir fq_dir

	gsutil ls ~{path_to_fq}
	gsutil -m cp -r ~{path_to_fq}* fq_dir/

	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --overlap_reads --meta --repo --adaptor_removal --merge --bbmerge_report
	Rscript /Code/BBMerge.R Report/Merge/ Report/

	cat Report/stdout.txt
	cat Report/stderr.txt

	ls Report/Merge/
	Rscript /Code/Contamination.R Report/Merge/ Report/ ~{path_to_flist} ~{joined_threshold} ~{contamination_threshold}
	tar -czvf Merge.tar.gz Results/Merge
	find . -type f
	#cat Results/stdout.txt
	#cat Results/stderr.txt
	>>>
	output {
		File rawfilelist = "Results/Fq_metadata/rawfilelist.tsv"
		File missing_files = "Results/missing_files.tsv" 
		File merge_tar = "Merge.tar.gz"
		File bbmergefields = "Report/Merge/bbmergefields.tsv"
		File BBmerge_performance_absolute_report = "Report/BBmerge_performance_absolute_report.svg"
		File BBmerge_performance_percentage_report = "Report/BBmerge_performance_percentage_report.svg"
		File BBmerge_performace_absolute_discarded = "Report/BBmerge_performace_absolute_discarded.svg"	
		File Barcode_report_abs = "Report/Barcode_report_abs.svg"
		File Barcode_report_per = "Report/Barcode_report_per.svg"
		File Insert_size = "Report/Insert_size.png"
		File Match_report_abs = "Report/Match_report_abs.svg"
		File Match_report_per = "Report/Match_report_per.svg"

		File barcodes_report_bbmerge = "Report/barcodes_report_bbmerge.tsv"
		File hamming_distances_forward = "Report/hamming_forward.tsv"
		File hamming_distances_reverse = "Report/hamming_reverse.tsv"	
		#?REMOVE
		#File config_MiSeq = config_json
		#File ampseq_bbmerge_process_stdout = stdout()
		#File ampseq_bbmerge_process_stderr = stderr()

	}
	runtime {
		cpu: 1
		memory: "1 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: 'jorgeamaya/ci_barcode_terra:v1'
	}
}
