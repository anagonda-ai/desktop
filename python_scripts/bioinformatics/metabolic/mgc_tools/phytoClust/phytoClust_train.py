import os
import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.optimizers import Adam

# Simulated dataset (Replace with real biological data)
X = np.array([
    [150, 1, 0.42, 8.5],
    [200, 0, 0.45, 7.2],
    [180, 1, 0.47, 6.8],
    [220, 1, 0.43, 8.9],
    [170, 0, 0.48, 7.5],
    [300, 1, 0.50, 9.0],
    [250, 0, 0.44, 6.5],
    [275, 1, 0.46, 8.3],
    [190, 0, 0.41, 7.9],
    [210, 1, 0.49, 9.2]
])

y_clusters = np.array([1, 0, 1, 1, 0, 1, 0, 1, 0, 1])  # 1 = Cluster, 0 = Not a Cluster
y_pathway = np.array([0, 1, 0, 2, 1, 0, 1, 2, 1, 0])  # Pathway labels

# Split into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y_clusters, test_size=0.2, random_state=42)
X_train_p, X_test_p, y_train_p, y_test_p = train_test_split(X, y_pathway, test_size=0.2, random_state=42)

# RandomForest Model for Cluster Detection
cluster_model = RandomForestClassifier(n_estimators=200, random_state=42)
cluster_model.fit(X_train, y_train)

# GradientBoosting Model for Pathway Prediction
pathway_model = GradientBoostingClassifier(n_estimators=200, learning_rate=0.05)
pathway_model.fit(X_train_p, y_train_p)

# Deep Learning Model for Pathway Prediction
deep_model = Sequential([
    Dense(128, activation='relu', input_shape=(X_train.shape[1],)),
    Dropout(0.3),
    Dense(64, activation='relu'),
    Dropout(0.2),
    Dense(len(set(y_pathway)), activation='softmax')
])

deep_model.compile(optimizer=Adam(learning_rate=0.001), loss='sparse_categorical_crossentropy', metrics=['accuracy'])
deep_model.fit(X_train_p, y_train_p, epochs=20, batch_size=4, validation_split=0.2, verbose=1)

current_file_path = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(current_file_path, "models")
os.makedirs(output_dir, exist_ok=True)

# Save models
joblib.dump(cluster_model, os.path.join(output_dir, "phytoClust_optimized_model.pkl"))
joblib.dump(pathway_model, os.path.join(output_dir, "pathway_prediction_model.pkl"))
deep_model.save(os.path.join(output_dir, "phytoClust_deep_model.h5"))

print("✅ Models trained and saved successfully!")
