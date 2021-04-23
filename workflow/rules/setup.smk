# An example collection of Snakemake rules imported in the main Snakefile.

rule install_software:
    # in case ./install_software.sh has not been run yet,
    # this rule can be used to do it automatically
    # this should only be used for binaries or software that can not
    # be downloaded through conda
    output:
        "bin/plink",
        "bin/plink2",
        "bin/gcta/gcta64",
        "bin/gctb/gctb"
    shell:
        "bash install_software.sh"


