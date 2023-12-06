alto cromwell run -s cumulus.aws.science.roche.com \
    -m cumulus/demultiplexing \
    -i ch10_souporcell_inputs.json \
    -o ch10_souporcell_inputs_updated.json \
    -b s3://gred-cumulus-data/AVITI_pilot/demultiplex_souporcell/upload \
    --profile gred-cumulus
