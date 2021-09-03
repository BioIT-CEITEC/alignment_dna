1. Velmi jednoducha dokumentacia (teda skoro ziadna) k execution na TES je tu: https://snakemake.readthedocs.io/en/stable/executing/cloud.html#executing-a-snakemake-workflow-via-ga4gh-tes
2. Aby snakemake so vsetkym co potrebuje fungoval, robila som full install podla: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html (vratane instalce mamby ale to asi mate)
3. Nasledne: `conda activate snakemake`
4. Treba vyexportovat tieto 2 premenne tak ako pise v dokumentacii: `export CONDA_ENVS_PATH=/tmp/conda` + `export CONDA_PKGS_DIRS=/tmp/conda`
5. Command vyzeral napriklad takto: `snakemake --tes https://tesk-na.cerit-sc.cz/ga4gh/tes --use-conda --envvars CONDA_PKGS_DIRS CONDA_ENVS_PATH --conda-prefix $CONDA_ENVS_PATH --jobs 8 -s Snakefile --verbose --reason --conda-frontend mamba`. Pridane su `--envvars`, url na TES instanciu. 


Note: `pip3 install py-tes` / `pip3 install tes` je potrebne vykonat predtym, nez sa pusti snakemake, ale v tom conda snakemake prostredi
Note2: Cesty k jednotlivym suborom som sa snazila co najviac generalizovat, pripadne nechat ako boli ale v niektorych miestach proste nesedeli a tak som ich napisala cele, napr. `input = { 'ref' : S3.remote(expand("cerit_sc_test/references/homsap/GRCh37-p13/index/BWA/GRCh37-p13.{ref}",ref=["amb","ann","bwt","pac","sa"]))}`
Note3: input a output by sa mali zjednotit tak, aby sa vzdy pouzivalo na konci indexovanie `[0]`, najlepsie rovno v rules aby to netrebalo pisat vsade do skriptov. Ja som to ale miesala podla toho, kde mi skakali errory. Dovod vid tu: https://github.com/snakemake/snakemake/issues/737
