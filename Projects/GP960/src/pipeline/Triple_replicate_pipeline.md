```bash
#!/bin/bash
## step1: pp_rawdata.py
# python ../data_processing/triple_replicate/01.pp_rawdata.py \
#     --rawdata_dir ../../data/raw/triple_replicate/NaCl \
#     --plate_file ../../data/raw/triple_replicate/plate-gene-mapping.xlsx \
#     --result_dir ../../data/results/triple_replicate/NaCl/01.ppraw_data \
#     --c1_label nonstress \
#     --c2_label stress

# ## step2: data_cleaning.py
# python ../data_processing/triple_replicate/02.data_cleaning.py \
#     --input_dir ../../data/results/triple_replicate/NaCl/01.ppraw_data \
#     --output_dir ../../data/results/triple_replicate/NaCl/02.cleaned_data \
#     --mapping_file ../../data/raw/triple_replicate/plate-gene-mapping.xlsx \
#     --threshold 0 \

# # ## step3: pp_racalculate_fitnesswdata.py
# python ../data_processing/triple_replicate/03.calculate_fitness.py \
#     --input_dir ../../data/results/triple_replicate/NaCl/02.cleaned_data \
#     --output_dir ../../data/results/triple_replicate/NaCl/03.stress_tolerance_analysis \
#     --rf_threshold_stress 1.05 \
#     --comprehensive_rf_threshold_stress 1.05 \
#     --rf_threshold_nonstress 0.95 \
#     --comprehensive_rf_threshold_nonstress 0.95 \

# ## step4: plot_gene_od.py
python ../data_processing/triple_replicate/04.plot_gene_od.py \
    --input_dir ../../data/results/triple_replicate/NaCl/02.cleaned_data \
    --output_dir ../../data/results/triple_replicate/NaCl/04.gene_growth_curves \
    --analysis_dir ../../data/results/triple_replicate/NaCl/03.stress_tolerance_analysis \

# python ../data_processing/triple_replicate/01.pp_rawdata.py --help
# python ../data_processing/triple_replicate/02.data_cleaning.py --help
# python ../data_processing/triple_replicate/03.calculate_fitness.py --help
# python ../data_processing/triple_replicate/04.plot_gene_od.py --help

```