
```bash
#!/bin/bash
## step1: pp_rawdata.py
# python ../data_processing/single_replicate/01.pp_rawdata.py \
#     --rawdata_dir ../../data/raw/single_replicate/highmethanol \
#     --plate_file ../../data/raw/single_replicate/plate-gene-mapping.xlsx \
#     --result_dir ../../data/results/single_replicate/highmethanol/01.ppraw_data \
#     --c1_label stress \
#     --c2_label nonstress

# ## step2: data_cleaning.py
# python ../data_processing/single_replicate/02.data_cleaning.py \
#     --input_dir ../../data/results/single_replicate/highmethanol/01.ppraw_data \
#     --output_dir ../../data/results/single_replicate/highmethanol/02.cleaned_data \
#     --mapping_file ../../data/raw/single_replicate/plate-gene-mapping.xlsx \
#     --threshold 1 \

# ## step3: pp_racalculate_fitnesswdata.py (optional)
# python ../data_processing/single_replicate/03.calculate_fitness.py \
#     --input_dir ../../data/results/single_replicate/highmethanol/02.cleaned_data \
#     --output_dir ../../data/results/single_replicate/highmethanol/03.stress_tolerance_analysis \
#     --rf_threshold_stress 1.05 \
#     --comprehensive_rf_threshold_stress 1.05 \
#     --rf_threshold_nonstress 0.95 \
#     --comprehensive_rf_threshold_nonstress 0.95 \

# ## step4: plot_gene_od.py (optional)
# python ../data_processing/single_replicate/04.plot_gene_od.py \
#     --input_dir ../../data/results/single_replicate/highmethanol/02.cleaned_data \
#     --output_dir ../../data/results/single_replicate/highmethanol/04.gene_growth_curves \
#     --analysis_dir ../../data/results/single_replicate/highmethanol/03.stress_tolerance_analysis \


## step5: calculate_full_features.py
# python ../data_processing/single_replicate/05.calculate_full_features.py \
#     --input_dir ../../data/results/single_replicate/highmethanol/02.cleaned_data \
#     --output_dir ../../data/results/single_replicate/highmethanol/05.full_features

## step6: plot_raw_correlation_heatmap_exact.py
# python ../data_processing/single_replicate/06.plot_raw_correlation_heatmap_exact.py \
#     --input ../../data/results/single_replicate/highmethanol/05.full_features/all_features_raw_50plus.csv \
#     --output_dir ../../data/results/single_replicate/highmethanol/06.plot_feature_selection \

## step6: plot_raw_feature_scatter.py
# python ../data_processing/single_replicate/06.plot_raw_feature_scatter.py \
#     --input ../../data/results/single_replicate/highmethanol/05.full_features/all_features_raw_50plus.csv \
#     --output_dir ../../data/results/single_replicate/highmethanol/06.plot_feature_selection \
#     --feat1 max_OD \
#     --feat2 AUC \

# step7: calculate_rf_matrix.py
# python ../data_processing/single_replicate/07.calculate_rf_matrix.py \
#     --input_features ../../data/results/single_replicate/highmethanol/05.full_features/all_features_raw_50plus.csv \
#     --input_selected ../../data/results/single_replicate/highmethanol/06.plot_feature_selection/selected_features_list.csv \
#     --output_dir ../../data/results/single_replicate/highmethanol/07.select_features_RF \

# step8: screen_top5_per_feature.py
# python ../data_processing/single_replicate/08.screen_top5_per_feature.py \
#     --input ../../data/results/single_replicate/highmethanol/07.select_features_RF/rf_matrix_final.csv \
#     --output_dir ../../data/results/single_replicate/highmethanol/08.top5_strains \
#     --top_n 5


# step9: plot_top_strains_curves.py
# python ../data_processing/single_replicate/09.plot_top_strains_curves.py \
#     --input_dir ../../data/results/single_replicate/highmethanol/02.cleaned_data \
#     --top5_file ../../data/results/single_replicate/highmethanol/08.top5_strains/top5_strains_per_feature.csv \
#     --output_dir ../../data/results/single_replicate/highmethanol/09.plot_top5_curves_by_feature \


# step10: 10.plot_dotplot.py
# python ../data_processing/single_replicate/10.plot_dotplot.py \
#     --top5_file ../../data/results/single_replicate/highmethanol/08.top5_strains/top5_strains_per_feature.csv \
#     --output_dir ../../data/results/single_replicate/highmethanol/10.dotplot \



python ../data_processing/single_replicate/01.pp_rawdata.py --help
python ../data_processing/single_replicate/02.data_cleaning.py --help
python ../data_processing/single_replicate/03.calculate_fitness.py --help ##(optional)
python ../data_processing/single_replicate/04.plot_gene_od.py --help ##(optional)
python ../data_processing/single_replicate/05.calculate_full_features.py --help
python ../data_processing/single_replicate/06.plot_raw_correlation_heatmap_exact.py --help
python ../data_processing/single_replicate/06.plot_raw_feature_scatter.py --help
python ../data_processing/single_replicate/07.calculate_rf_matrix.py --help
python ../data_processing/single_replicate/08.screen_top5_per_feature.py --help
python  ../data_processing/single_replicate/09.plot_top_strains_curves.py --help
python ../data_processing/single_replicate/10.plot_dotplot.py --help

```