import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.preprocessing import StandardScaler
import seaborn as sns

# Set plotting style
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ==============================================================================
# 1. Load data
# ==============================================================================

df_135 = pd.read_csv('new_135.csv', index_col=0)
df_297 = pd.read_csv('new_297.csv', index_col=0)
df_658 = pd.read_csv('new_658.csv', index_col=0)

print("=" * 60)
print("Data Overview")
print("=" * 60)
print(f"135 training set: {df_135.shape}")
print(f"297 test set: {df_297.shape}")
print(f"658 test set: {df_658.shape}")
print()

print("135 training set columns:", df_135.columns.tolist())
print("135 training set head:")
print(df_135.head())
print()

print("Type column value distribution (135 training set):")
print(df_135['Type'].value_counts())
print()

print("Data statistics:")
print(df_135[['AS', 'EXPORT', 'TE', 'DT']].describe())

# ==============================================================================
# 2. Data visualization
# ==============================================================================

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

features = ['AS', 'EXPORT', 'TE']
target = 'DT'

# Feature distributions
for i, feat in enumerate(features):
    axes[0, i].hist(df_135[feat].dropna(), bins=20, edgecolor='black', alpha=0.7)
    axes[0, i].set_title(f'{feat} distribution (n={df_135[feat].notna().sum()})')
    axes[0, i].set_xlabel(feat)
    axes[0, i].set_ylabel('Frequency')

# Target distribution
axes[1, 0].hist(df_135[target].dropna(), bins=20, edgecolor='black', alpha=0.7, color='green')
axes[1, 0].set_title(f'{target} distribution')
axes[1, 0].set_xlabel(target)
axes[1, 0].set_ylabel('Frequency')

# Correlation heatmap
corr_data = df_135[features + [target]].corr()
sns.heatmap(corr_data, annot=True, cmap='coolwarm', center=0, ax=axes[1, 1])
axes[1, 1].set_title('Feature-target correlations')

# AS vs DT scatter plot colored by Type
types = df_135['Type'].unique()
colors = plt.cm.Set1(np.linspace(0, 1, len(types)))
for t, color in zip(types, colors):
    subset = df_135[df_135['Type'] == t]
    axes[1, 2].scatter(subset['AS'], subset['DT'], alpha=0.6, c=[color], label=t)
axes[1, 2].set_xlabel('AS (PSI)')
axes[1, 2].set_ylabel('DT')
axes[1, 2].set_title('AS vs DT (colored by splicing type)')
axes[1, 2].legend()

plt.tight_layout()
plt.savefig('data_distribution.png', dpi=150, bbox_inches='tight')
plt.show()

print("\nCorrelation matrix:")
print(corr_data)

# ==============================================================================
# 3. AS distribution by Type (supplementary analysis)
# ==============================================================================

print("\n" + "=" * 60)
print("AS distribution by splicing type")
print("=" * 60)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for t in types:
    subset = df_135[df_135['Type'] == t]
    axes[0].hist(subset['AS'], bins=15, alpha=0.5, label=t, edgecolor='black')
axes[0].set_title('AS distribution by splicing type')
axes[0].set_xlabel('AS (PSI)')
axes[0].set_ylabel('Frequency')
axes[0].legend()

df_135.boxplot(column='AS', by='Type', ax=axes[1])
axes[1].set_title('AS boxplot by splicing type')
axes[1].set_xlabel('Splicing type')
axes[1].set_ylabel('AS (PSI)')

plt.tight_layout()
plt.savefig('as_by_type.png', dpi=150, bbox_inches='tight')
plt.show()

print("\nAS statistics by splicing type:")
as_stats = df_135.groupby('Type')['AS'].agg(['count', 'mean', 'std', 'min', 'max'])
print(as_stats)

# ==============================================================================
# 4. Data preprocessing
# ==============================================================================

print("\n" + "=" * 60)
print("Data preprocessing")
print("=" * 60)

print("\nMissing value check:")
print(f"135 set missing values:\n{df_135[features + [target]].isnull().sum()}")

df_135_clean = df_135[features + [target]].dropna()
df_297_clean = df_297[features + [target]].dropna()
df_658_clean = df_658[features + [target]].dropna()

print(f"\nCleaned 135 set: {df_135_clean.shape}")
print(f"Cleaned 297 set: {df_297_clean.shape}")
print(f"Cleaned 658 set: {df_658_clean.shape}")

X_train = df_135_clean[features].values
y_train = df_135_clean[target].values

X_297 = df_297_clean[features].values
y_297 = df_297_clean[target].values

X_658 = df_658_clean[features].values
y_658 = df_658_clean[target].values

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_297_scaled = scaler.transform(X_297)
X_658_scaled = scaler.transform(X_658)

# ==============================================================================
# 5. Model training - Linear Regression
# ==============================================================================

print("\n" + "=" * 60)
print("Model training - Linear Regression")
print("=" * 60)

# 5.1 Ordinary Linear Regression
lr = LinearRegression()
lr.fit(X_train, y_train)

print("\n[Ordinary Linear Regression - original scale]")
print(f"AS weight: {lr.coef_[0]:.4f}")
print(f"EXPORT weight: {lr.coef_[1]:.4f}")
print(f"TE weight: {lr.coef_[2]:.4f}")
print(f"Intercept: {lr.intercept_:.4f}")
print(f"\nFormula: DT = {lr.coef_[0]:.4f} * AS + {lr.coef_[1]:.4f} * EXPORT + {lr.coef_[2]:.4f} * TE + {lr.intercept_:.4f}")

cv_scores = cross_val_score(lr, X_train, y_train, cv=5, scoring='r2')
print(f"\n5-fold CV R2 mean: {cv_scores.mean():.4f} +/- {cv_scores.std():.4f}")

y_train_pred = lr.predict(X_train)
train_r2 = r2_score(y_train, y_train_pred)
train_mse = mean_squared_error(y_train, y_train_pred)
train_mae = mean_absolute_error(y_train, y_train_pred)

print(f"\nTraining set performance:")
print(f"  R2: {train_r2:.4f}")
print(f"  MSE: {train_mse:.4f}")
print(f"  MAE: {train_mae:.4f}")

# 5.2 Standardized Linear Regression
lr_scaled = LinearRegression()
lr_scaled.fit(X_train_scaled, y_train)

print("\n[Standardized Linear Regression - comparable weights]")
print(f"AS standardized weight: {lr_scaled.coef_[0]:.4f}")
print(f"EXPORT standardized weight: {lr_scaled.coef_[1]:.4f}")
print(f"TE standardized weight: {lr_scaled.coef_[2]:.4f}")

# 5.3 Ridge Regression
ridge = Ridge(alpha=1.0)
ridge.fit(X_train_scaled, y_train)
print(f"\n[Ridge Regression (alpha=1) - standardized weights]")
print(f"AS: {ridge.coef_[0]:.4f}, EXPORT: {ridge.coef_[1]:.4f}, TE: {ridge.coef_[2]:.4f}")

# 5.4 Lasso Regression
lasso = Lasso(alpha=0.01)
lasso.fit(X_train_scaled, y_train)
print(f"\n[Lasso Regression (alpha=0.01) - standardized weights]")
print(f"AS: {lasso.coef_[0]:.4f}, EXPORT: {lasso.coef_[1]:.4f}, TE: {lasso.coef_[2]:.4f}")

# ==============================================================================
# 6. Test set evaluation
# ==============================================================================

print("\n" + "=" * 60)
print("Test set evaluation")
print("=" * 60)

y_297_pred = lr.predict(X_297)
y_658_pred = lr.predict(X_658)

r2_297 = r2_score(y_297, y_297_pred)
r2_658 = r2_score(y_658, y_658_pred)

mse_297 = mean_squared_error(y_297, y_297_pred)
mse_658 = mean_squared_error(y_658, y_658_pred)

mae_297 = mean_absolute_error(y_297, y_297_pred)
mae_658 = mean_absolute_error(y_658, y_658_pred)

print(f"\n[297 gene set]")
print(f"  Samples: {len(y_297)}")
print(f"  R2: {r2_297:.4f}")
print(f"  MSE: {mse_297:.4f}")
print(f"  MAE: {mae_297:.4f}")

print(f"\n[658 gene set]")
print(f"  Samples: {len(y_658)}")
print(f"  R2: {r2_658:.4f}")
print(f"  MSE: {mse_658:.4f}")
print(f"  MAE: {mae_658:.4f}")

print(f"\n[Comparison]")
print(f"  297 R2 / Training R2 = {r2_297 / train_r2:.4f}")
print(f"  658 R2 / Training R2 = {r2_658 / train_r2:.4f}")

# ==============================================================================
# 7. Predicted vs actual plots
# ==============================================================================

type_train = df_135.loc[df_135_clean.index, 'Type'] if 'Type' in df_135.columns else None
type_297 = df_297.loc[df_297_clean.index, 'Type'] if 'Type' in df_297.columns else None
type_658 = df_658.loc[df_658_clean.index, 'Type'] if 'Type' in df_658.columns else None

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

datasets = [
    (y_train, y_train_pred, type_train, '135 training set', train_r2),
    (y_297, y_297_pred, type_297, '297 test set', r2_297),
    (y_658, y_658_pred, type_658, '658 test set', r2_658)
]

for ax, (y_true, y_pred, types, title, r2) in zip(axes, datasets):
    if types is not None:
        unique_types = types.unique()
        colors = plt.cm.Set1(np.linspace(0, 1, len(unique_types)))
        for t, color in zip(unique_types, colors):
            mask = types == t
            if mask.sum() > 0:
                ax.scatter(y_true[mask], y_pred[mask], alpha=0.6, c=[color],
                          label=t, edgecolors='k', linewidth=0.5)
        ax.legend(loc='best', fontsize=8)
    else:
        ax.scatter(y_true, y_pred, alpha=0.5, c='blue', edgecolors='k', linewidth=0.5)

    min_val = min(y_true.min(), y_pred.min())
    max_val = max(y_true.max(), y_pred.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='y=x')

    ax.set_xlabel('Actual')
    ax.set_ylabel('Predicted')
    ax.set_title(f'{title}\nR2 = {r2:.4f}')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('prediction_comparison.png', dpi=150, bbox_inches='tight')
plt.show()

# ==============================================================================
# 8. Residual analysis
# ==============================================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

residuals_train = y_train - y_train_pred
residuals_297 = y_297 - y_297_pred
residuals_658 = y_658 - y_658_pred

axes[0].hist(residuals_train, bins=20, edgecolor='black', alpha=0.7, color='blue')
axes[0].axvline(x=0, color='red', linestyle='--')
axes[0].set_title(f'135 set residuals\nmean={residuals_train.mean():.4f}')
axes[0].set_xlabel('Residual')

axes[1].hist(residuals_297, bins=20, edgecolor='black', alpha=0.7, color='orange')
axes[1].axvline(x=0, color='red', linestyle='--')
axes[1].set_title(f'297 set residuals\nmean={residuals_297.mean():.4f}')
axes[1].set_xlabel('Residual')

axes[2].hist(residuals_658, bins=20, edgecolor='black', alpha=0.7, color='red')
axes[2].axvline(x=0, color='red', linestyle='--')
axes[2].set_title(f'658 set residuals\nmean={residuals_658.mean():.4f}')
axes[2].set_xlabel('Residual')

plt.tight_layout()
plt.savefig('residual_analysis.png', dpi=150, bbox_inches='tight')
plt.show()

# ==============================================================================
# 9. Feature weights visualization
# ==============================================================================

fig, ax = plt.subplots(figsize=(8, 5))

features_name = ['AS', 'EXPORT', 'TE']
weights = lr.coef_
colors = ['#2ecc71', '#3498db', '#e74c3c']

bars = ax.bar(features_name, weights, color=colors, edgecolor='black')
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax.set_ylabel('Weight coefficient')
ax.set_title('Feature weights for DT prediction')
ax.grid(axis='y', alpha=0.3)

for bar, w in zip(bars, weights):
    height = bar.get_height()
    ax.annotate(f'{w:.4f}',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3 if height >= 0 else -15),
                textcoords="offset points",
                ha='center', va='bottom')

plt.tight_layout()
plt.savefig('feature_weights.png', dpi=150, bbox_inches='tight')
plt.show()

# ==============================================================================
# 10. Type effect on prediction error
# ==============================================================================

print("\n" + "=" * 60)
print("Splicing type effect on prediction error")
print("=" * 60)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

if type_train is not None:
    residuals_train_abs = np.abs(residuals_train)
    train_df = pd.DataFrame({'Type': type_train, 'Abs_Residual': residuals_train_abs})

    train_df.boxplot(column='Abs_Residual', by='Type', ax=axes[0])
    axes[0].set_title('Training set: absolute residual by splicing type')
    axes[0].set_xlabel('Splicing type')
    axes[0].set_ylabel('Absolute residual')

    print("\nTraining set - prediction error by splicing type:")
    error_stats = train_df.groupby('Type')['Abs_Residual'].agg(['count', 'mean', 'std'])
    print(error_stats)

if type_297 is not None:
    residuals_297_abs = np.abs(residuals_297)
    test297_df = pd.DataFrame({'Type': type_297, 'Abs_Residual': residuals_297_abs})

    test297_df.boxplot(column='Abs_Residual', by='Type', ax=axes[1])
    axes[1].set_title('297 test set: absolute residual by splicing type')
    axes[1].set_xlabel('Splicing type')
    axes[1].set_ylabel('Absolute residual')

    print("\n297 test set - prediction error by splicing type:")
    error_stats_297 = test297_df.groupby('Type')['Abs_Residual'].agg(['count', 'mean', 'std'])
    print(error_stats_297)

plt.tight_layout()
plt.savefig('error_by_type.png', dpi=150, bbox_inches='tight')
plt.show()

# ==============================================================================
# 11. Summary report
# ==============================================================================

print("\n" + "=" * 60)
print("Summary report")
print("=" * 60)

print(f"""
1. Model formula (based on 135 training genes):
   DT = {lr.coef_[0]:.4f} * AS + {lr.coef_[1]:.4f} * EXPORT + {lr.coef_[2]:.4f} * TE + {lr.intercept_:.4f}

2. Feature importance (absolute weight):
   - Translation efficiency (TE): {abs(lr.coef_[2]):.4f}
   - Nuclear export (EXPORT): {abs(lr.coef_[1]):.4f}
   - Splicing (AS): {abs(lr.coef_[0]):.4f}

   -> {'TE' if abs(lr.coef_[2]) > abs(lr.coef_[1]) and abs(lr.coef_[2]) > abs(lr.coef_[0]) else 'EXPORT' if abs(lr.coef_[1]) > abs(lr.coef_[0]) else 'AS'} has the strongest effect on DT

3. Model performance:
   - 135 training set R2: {train_r2:.4f}
   - 297 test set R2: {r2_297:.4f}
   - 658 test set R2: {r2_658:.4f}

   -> Model fits {'297' if r2_297 > r2_658 else '658'} test set better
   -> 297 R2 / training R2 ratio: {r2_297 / train_r2:.4f}
   -> 658 R2 / training R2 ratio: {r2_658 / train_r2:.4f}

4. Type column (splicing event type) was used for visualization only:
   - AS distribution by type shown in as_by_type.png
   - Prediction error by type shown in error_by_type.png
""")

results = {
    'AS_weight': lr.coef_[0],
    'EXPORT_weight': lr.coef_[1],
    'TE_weight': lr.coef_[2],
    'intercept': lr.intercept_,
    'train_r2': train_r2,
    'r2_297': r2_297,
    'r2_658': r2_658,
    'mse_297': mse_297,
    'mse_658': mse_658
}

pd.Series(results).to_csv('model_results.csv')
print("\nModel results saved to model_results.csv")
