dx generate_batch_inputs -istage-G0q44PQ433GYV97qKk0zkb4J.fastqs='(.*)_L(.*)_R1(.*)' -istage-G0q44PQ433GYV97qKk0zkb4J.fastqs2='(.*)_L(.*)_R2(.*)' -istage-G0qpXy0433Gv75XbPJ3xj8jV.reads_fastqgzs='(.*)_L(.*)_R1(.*)' -istage-G0qpXy0433Gv75XbPJ3xj8jV.reads2_fastqgzs='(.*)_L(.*)_R2(.*)'


head -n 1 dx_batch.0000.tsv > temp.tsv && tail -n +2 dx_batch.0000.tsv | awk '{sub($6, "["$6, $6); sub($7, $7"]", $7); sub($8, "["$8"]", $8); sub($9, "["$9"]", $9); print}' >> temp.tsv; tr -d '\r' < temp.tsv > b1.tsv; rm temp.tsv

dx run project-Fkb6Gkj433GVVvj73J7x8KbV:workflow-G18JB78433Gbf8YJ9BBF2PQZ --batch-tsv b2.tsv --destination "project-G18J7j04fy6jgz3937gJ9JfG:/output/uranus_210316_1"

dx run project-Fkb6Gkj433GVVvj73J7x8KbV:applet-Fz93FfQ433Gvf6pKFZYbXZQf -ieggd_multiqc_config_file="project-Fkb6Gkj433GVVvj73J7x8KbV:file-G0K191j433Gv6JG63b43z8Gy" -iproject_for_multiqc="003_210312_K00178_0315_AHL5GMBBXY2" -iss_for_multiqc="uranus_210316_1" -icustom_coverage="true" --destination "project-G18J7j04fy6jgz3937gJ9JfG:/output/uranus_210316_1/MultiQC_v1.1.1"


dx generate_batch_inputs -istage-G0QQ8jj433Gxyx2K8xfPyV7B.input_vcf='(.*).vcf.gz$' --path=project-G18J7j04fy6jgz3937gJ9JfG:/output/uranus_210316_1/pindel_filtering_v1.0.0
# Run annotation workflow for pindel
dx run project-Fkb6Gkj433GVVvj73J7x8KbV:workflow-G0bj928433GxbZY9632xXP9J --batch-tsv dx_batch.0000.tsv --destination "project-G18J7j04fy6jgz3937gJ9JfG:/output/uranus_210316_1/pindel_annotation"


dx generate_batch_inputs -istage-G0QQ8jj433Gxyx2K8xfPyV7B.input_vcf='(.*).vcf.gz$' --path=project-G18J7j04fy6jgz3937gJ9JfG:/output/uranus_210316_1/cgppindel_filtering_v1.0.0
# Run annotation workflow for cgppindel
dx run project-Fkb6Gkj433GVVvj73J7x8KbV:workflow-G0bj928433GxbZY9632xXP9J --batch-tsv dx_batch.0000.tsv --destination "project-G18J7j04fy6jgz3937gJ9JfG:/output/uranus_210316_1/cgppindel_annotation"