# Models Directory

This directory should contain the trained machine learning model files.

## Required Files

- `combined_model.pkl` - Trained XGBoost model (~258 KB)
- `combined_model_feature_order.txt` - Feature order file (required for predictions)
- `combined_all_features_weights.csv` - Feature importances (optional, for reference)

## Obtaining the Files

The model files should be obtained from:
- Original location: `/groups/itay_mayrose/alongonda/desktop/mibig_validate/combined/`

To copy the files:

```bash
cp /groups/itay_mayrose/alongonda/desktop/mibig_validate/combined/combined_model.pkl FINAL_TOOL/models/
cp /groups/itay_mayrose/alongonda/desktop/mibig_validate/combined/combined_model_feature_order.txt FINAL_TOOL/models/
cp /groups/itay_mayrose/alongonda/desktop/mibig_validate/combined/combined_all_features_weights.csv FINAL_TOOL/models/
```

Alternatively, run the setup script which will help you verify and prepare all required files:

```bash
python FINAL_TOOL/setup.py
```

