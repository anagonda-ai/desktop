import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import re

# נתיב לתיקייה
directory_path = "/groups/itay_mayrose_nosnap/alongonda/Plant_MGC/unique_clusters_sliding_window_outputs"

# רשימה לאחסון window_size ומספר השורות
window_sizes = []
row_counts = []

# מעבר על כל הקבצים בתיקייה
for file_name in os.listdir(directory_path):
    file_path = os.path.join(directory_path, file_name)
    if os.path.isfile(file_path):
        # חילוץ הספרה של ה-window_size מתוך שם הקובץ
        match = re.search(r'window_size_(\d+)', file_name)
        if match:
            window_size = int(match.group(1))  # המרת הספרה למספר שלם
            try:
                # קריאת הקובץ
                df = pd.read_csv(file_path, sep="\t")  # הנחה שהקבצים הם Tab-delimited
                window_sizes.append(window_size)
                row_counts.append(len(df))
            except Exception as e:
                print(f"Could not process file {file_name}: {e}")

# יצירת היסטוגרמה
plt.figure(figsize=(10, 6))
bars = plt.bar(window_sizes, row_counts, color='skyblue')
plt.xticks(window_sizes)
plt.xlabel('Window Size')
plt.ylabel('Number of Candidates (Rows)')
plt.title('Distribution of Candidates Across Window Sizes')
plt.tight_layout()

# הוספת ערכים ספציפיים לכל עמודה
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval, int(yval), va='bottom')  # va='bottom' to place text above the bar

# שמירת ההיסטוגרמה לקובץ
output_path = "/groups/itay_mayrose_nosnap/alongonda/Plant_MGC/unique_clusters_sliding_window_outputs/window_size_histogram.png"
plt.savefig(output_path, dpi=300)  # שמירה ברזולוציה גבוהה
plt.close()

print(f"Histogram saved to {output_path}")
print(f"Current working directory: {os.getcwd()}")
print(f"Files in current directory: {os.listdir('.')}")

# פונקציה להצגת הקורלציה בין גודל החלון למספר המועמדים
def show_correlation(window_sizes, row_counts):
    correlation, _ = pearsonr(window_sizes, row_counts)
    print(f"Correlation between window size and number of candidates: {correlation:.2f}")

# הצגת הקורלציה
show_correlation(window_sizes, row_counts)