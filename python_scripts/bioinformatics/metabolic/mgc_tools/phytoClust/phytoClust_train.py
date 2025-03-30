import os
import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.optimizers import Adam

# Load real features (4 columns: length, is_core, gc_content, expression)
X = np.load("/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/metabolic/mgc_tools/phytoClust/arabidopsis_expression_mean_tpm.npy")

# Temporary dummy labels (you should replace this with real cluster/pathway labels)
y_clusters = np.random.randint(0, 2, X.shape[0])
y_pathway = np.random.randint(0, 3, X.shape[0])

# Split for both tasks
X_train, X_test, y_train, y_test = train_test_split(X, y_clusters, test_size=0.2, random_state=42)
X_train_p, X_test_p, y_train_p, y_test_p = train_test_split(X, y_pathway, test_size=0.2, random_state=42)

# Model 1: Cluster classifier
cluster_model = RandomForestClassifier(n_estimators=200, random_state=42)
cluster_model.fit(X_train, y_train)

# Model 2: Pathway classifier
pathway_model = GradientBoostingClassifier(n_estimators=200, learning_rate=0.05)
pathway_model.fit(X_train_p, y_train_p)

# Model 3: Deep Learning pathway predictor
deep_model = Sequential([
    Dense(128, activation='relu', input_shape=(X.shape[1],)),
    Dropout(0.3),
    Dense(64, activation='relu'),
    Dropout(0.2),
    Dense(len(set(y_pathway)), activation='softmax')
])
deep_model.compile(optimizer=Adam(0.001), loss='sparse_categorical_crossentropy', metrics=['accuracy'])
deep_model.fit(X_train_p, y_train_p, epochs=20, batch_size=8, validation_split=0.2, verbose=1)

# Save models
model_dir = os.path.join(os.path.dirname(__file__), "models")
os.makedirs(model_dir, exist_ok=True)
joblib.dump(cluster_model, os.path.join(model_dir, "phytoClust_optimized_model.pkl"))
joblib.dump(pathway_model, os.path.join(model_dir, "pathway_prediction_model.pkl"))
deep_model.save(os.path.join(model_dir, "phytoClust_deep_model.h5"))

print("✅ Models trained and saved successfully!")
