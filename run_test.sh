nextflow run \
        -resume \
        --input ./assets/test_samples.csv \
        --output ./data/output/test \
        --grouper_config ./assets/grouper_config.yaml \
        main.nf
