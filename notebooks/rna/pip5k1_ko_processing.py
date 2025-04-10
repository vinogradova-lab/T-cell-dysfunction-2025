import os

import glob

import multiprocessing


def download_gencode_hg38_and_build_index_star():
    os.system("sbatch generate_genome_index.sh")


def write_and_launch_bash_script(
    commands, working_dir, job_name, bash_file_directory, cores, partition
):
    text = [
        "#!/bin/bash",
        "#SBATCH --job-name={}".format(job_name),
        "#SBATCH --error={}/logs/{}.err".format(working_dir, job_name),
        "#SBATCH --output={}/logs/{}.out".format(working_dir, job_name),
        "#SBATCH -n {}".format(cores),
        "#SBATCH -p {}".format(partition),
        commands,
    ]

    open_script = open("{}/{}.sh".format(bash_file_directory, job_name), "w+")

    open_script.write("\n".join(text))

    open_script.close()

    output_id = (
        os.popen("sbatch < {}/{}.sh".format(bash_file_directory, job_name))
        .read()
        .split(" ")[-1]
        .strip()
    )

    return output_id


def run_processing_and_alignment(working_dir, input_directory, index_dir):
    for i in [
        "star_output",
        "filtered_star_output",
        "sorted_bam_files",
        "bash_scripts",
        "logs",
        "trimmed_fastqs",
    ]:

        if not os.path.exists("{}".format(i)):

            os.system("mkdir {}".format(i))

    for i in glob.glob("{}/*_R1_001.fastq.gz".format(input_directory)):

        print(i)

        file_prefix = i.split("/")[-1].strip().split("_R1_001.fastq.gz")[0]

        print(file_prefix)

        file_1 = "{}/{}_R1_001.fastq.gz".format(input_directory, file_prefix)

        print(file_1)

        file_1_prefix = file_1.split("/")[-1].split(".fastq")[0]

        file_2 = "{}/{}_R2_001.fastq.gz".format(input_directory, file_prefix)

        print(file_2)

        file_2_prefix = file_2.split("/")[-1].split(".fastq")[0]

        commands = [
            "cd {}".format(working_dir),
            "source /ru-auth/local/home/hsanford/scratch/miniconda3/etc/profile.d/conda.sh",
            "conda activate sample_align",
            "trim_galore --quality 15 -j 6 --fastqc --gzip --paired {} {} -o trimmed_fastqs".format(
                file_1, file_2
            ),
            "STAR --runMode alignReads --runThreadN 8 --genomeDir {} --twopassMode Basic --readFilesIn trimmed_fastqs/{}_val_1.fq.gz trimmed_fastqs/{}_val_2.fq.gz --outFileNamePrefix star_output/{} --readFilesCommand gunzip -c".format(
                index_dir, file_1_prefix, file_2_prefix, file_prefix
            ),
            "samtools view -@ 8 -bo filtered_star_output/{}.bam star_output/{}Aligned.out.sam".format(
                file_prefix, file_prefix
            ),
            "samtools sort filtered_star_output/{}.bam -o sorted_bam_files/{}_sorted.bam".format(
                file_prefix, file_prefix
            ),
            "samtools index -b sorted_bam_files/{}_sorted.bam sorted_bam_files/{}_sorted.bam.bai".format(
                file_prefix, file_prefix
            ),
        ]

        write_and_launch_bash_script(
            job_name = "pipeline_run_{}".format(file_prefix),
            commands = "\n".join(commands),
            bash_file_directory="bash_scripts",
            working_dir=working_dir,
            cores=8,
            partition="hpc_a10",
        )


class featureCountsAnalysis:

    def __init__(self, sample_annotation_dictionary, output_prefix, gtf_file):
        for i in [
            "results_deseq2",
            "pca_plots_deseq2",
            "dot_plots_deseq2",
            "gsea_plots_fgsea",
            "gsea_tables_fgsea",
        ]:

            if not os.path.exists(i):

                os.system("mkdir {}".format(i))

        self.working_dir = os.getcwd()

        self.output_prefix = output_prefix

        self.sample_annotation_dictionary = sample_annotation_dictionary

        self.file_prefixes = sample_annotation_dictionary.keys()

        self.matrix_col_info_path = "{}_header_info.txt".format(output_prefix)

        self.gtf_file = gtf_file

        self.matrix_path = "{}_matrix.txt".format(output_prefix)

    def generate_metainfo(self):

        # generate col_info for deseq2
        meta_data_lines = ["Sample\tCondition"]

        for i in self.file_prefixes:
            meta_data_lines.append(
                "{}\t{}".format(i, self.sample_annotation_dictionary[i])
            )

        open_matrix_col_info = open(self.matrix_col_info_path, "w+")

        open_matrix_col_info.write("\n".join(meta_data_lines))

        open_matrix_col_info.write("\n")

        open_matrix_col_info.close()

    def generate_feature_counts_matrix(self):
        commands = [
            "module load conda/4.9.2 && source activate featurecounts && featureCounts -p  -T 8 -F GTF -t exon -g gene_name -a {} -o {} sorted_bam_files/*bam".format(
                self.gtf_file, self.matrix_path
            ),
            'grep -v "#" {} > temp.txt && mv temp.txt {}'.format(
                self.matrix_path, self.matrix_path
            ),
        ]

        for i in commands:
            os.system(i)

    def generate_pca(self):
        pca_path = "pca_plots_deseq2/{}_pca.pdf".format(self.output_prefix)

        os.system(
            "module load R/R-4.1.2 && Rscript generate_pca.R {} {} {} {}".format(
                self.matrix_path,
                self.matrix_col_info_path,
                self.output_prefix,
                pca_path,
            )
        )

    def generate_contrast_deseq2(self, numerator, denominator):

        comparison_prefix = "{}_{}_{}".format(
            self.output_prefix, numerator, denominator
        )

        results_path = "results_deseq2/{}_results.txt".format(comparison_prefix)

        os.system(
            "module load R/R-4.1.2 && Rscript deseq_differential_analysis.R {} {} {} {} {}".format(
                self.matrix_path,
                self.matrix_col_info_path,
                numerator,
                denominator,
                results_path,
            )
        )

        dot_plot_path = "dot_plots_deseq2/{}_dot_plot.pdf".format(comparison_prefix)

        os.system(
            "module load R/R-4.1.2 && Rscript dot_plots.R {} {} {}".format(
                results_path, comparison_prefix, dot_plot_path
            )
        )

    def generate_gsea(self, deseq_table, species, gmt_file):
        name = deseq_table.split("/")[-1].split(".")[0]

        pdf_name = "gsea_plots_fgsea/{}.pdf".format(name)

        txt_file_name = "gsea_tables_fgsea/{}.txt".format(name)

        title = "GSEA in {} across {}".format(name, gmt_file).replace(" ", "_")

        os.system(
            "module load R/R-4.1.2 && Rscript gsea.R {} {} {} {} {} {}".format(
                deseq_table, species, gmt_file, title, pdf_name, txt_file_name
            )
        )


def main():

    # download_gencode_hg38_and_build_index_star()

    # data_path = "/igo/delivery/share/vardhans/Project_13524/MICHELLE_0548"
    # merge_lanes(data_path, "formatted_fastqs")
    run_processing_and_alignment(os.getcwd(), "pip5k1ko", "/ru-auth/local/home/hsanford/scratch/shared_reference/gencode_v41/gencode_v41_index")

    metadata = {

"Donor1-1-D8A-Ctrl":"D8A-Ctrl",
"Donor1-2-D8A-1A":"D8A-1A",
"Donor1-3-D8A-1B":"D8A-1B",
"Donor1-4-D8A-1C":"D8A-1C",
"Donor1-5-D8C-Ctrl":"D8C-Ctrl",
"Donor1-6-D8C-1A":"D8C-1A",
"Donor1-7-D8C-1B":"D8C-1B",
"Donor1-8-D8C-1C":"D8C-1C",
"Donor1-9-D8A-Ctrl":"D8A-Ctrl",
"Donor1-10-D8A-1A":"D8A-1A",
"Donor1-11-D8A-1B":"D8A-1B",
"Donor1-12-D8A-1C":"D8A-1C",
"Donor1-13-D8C-Ctrl":"D8C-Ctrl",
"Donor1-14-D8C-1A":"D8C-1A",
"Donor1-15-D8C-1B":"D8C-1B",
"Donor1-16-D8C-1C":"D8C-1C",
"Donor2-1-D8A-Ctrl":"D8A-Ctrl",
"Donor2-2-D8A-1A":"D8A-1A",
"Donor2-3-D8A-1B":"D8A-1B",
"Donor2-4-D8A-1C":"D8A-1C",
"Donor2-5-D8C-Ctrl":"D8C-Ctrl",
"Donor2-6-D8C-1A":"D8C-1A",
"Donor2-7-D8C-1B":"D8C-1B",
"Donor2-8-D8C-1C":"D8C-1C",
"Donor2-9-D8A-Ctrl":"D8A-Ctrl",
"Donor2-10-D8A-1A":"D8A-1A",
"Donor2-11-D8A-1B":"D8A-1B",
"Donor2-12-D8A-1C":"D8A-1C",
"Donor2-13-D8C-Ctrl":"D8C-Ctrl",
"Donor2-14-D8C-1A":"D8C-1A",
"Donor2-15-D8C-1B":"D8C-1B",
"Donor2-16-D8C-1C":"D8C-1C",
"Donor3-1-D8A-Ctrl":"D8A-Ctrl",
"Donor3-2-D8A-1A":"D8A-1A",
"Donor3-3-D8A-1B":"D8A-1B",
"Donor3-4-D8A-1C":"D8A-1C",
"Donor3-5-D8C-Ctrl":"D8C-Ctrl",
"Donor3-6-D8C-1A":"D8C-1A",
"Donor3-7-D8C-1B":"D8C-1B",
"Donor3-8-D8C-1C":"D8C-1C",
"Donor3-9-D8A-Ctrl":"D8A-Ctrl",
"Donor3-10-D8A-1A":"D8A-1A",
"Donor3-11-D8A-1B":"D8A-1B",
"Donor3-12-D8A-1C":"D8A-1C",
"Donor3-13-D8C-Ctrl":"D8C-Ctrl",
"Donor3-14-D8C-1A":"D8C-1A",
"Donor3-15-D8C-1B":"D8C-1B",
"Donor3-16-D8C-1C":"D8C-1C",
}
    
    pip5k1ko_analysis = featureCountsAnalysis(sample_annotation_dictionary=metadata,
                                              output_prefix="pip5k1ko", gtf_file="gencode.v41.chr_patch_hapl_scaff.annotation.gtf")


    pip5k1ko_analysis.generate_metainfo()
    pip5k1ko_analysis.generate_feature_counts_matrix()
    pip5k1ko_analysis.generate_pca()

    # p13524_analysis.generate_contrast_deseq2("D4-Ac", "D2")
    # p13524_analysis.generate_contrast_deseq2("D4-Chr", "D2")
    # p13524_analysis.generate_contrast_deseq2("D8-Ac", "D2")
    # p13524_analysis.generate_contrast_deseq2("D8-Chr", "D2")

    # for i in glob.glob("results_deseq2/*"):
    # 	p13524_analysis.generate_gsea(i, "Human", "h.all.v7.5.1.symbols.gmt")


if __name__ == "__main__":
    main()


rm *_R1_001*