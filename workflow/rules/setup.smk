# An example collection of Snakemake rules imported in the main Snakefile.

rule install_software:
    output:
        "bin/plink",
        "bin/plink2",
        "bin/gcta/gcta64",
        "bin/gctb/gctb"
    shell:
        "bash install_software.sh"
