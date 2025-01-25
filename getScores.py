import re
from pathlib import Path
import pandas as pd

def extract_as_score(sam_file):
    with open(sam_file, 'r') as f:
        for line in f:
            if match := re.search(r'AS:i:(-?\d+)', line):
                return int(match.group(1))
    return None

def extract_k(filepath):
    # Extract k value from path like alignment_results/reference/k_value/...
    k = filepath.parts[-2]  # Get the parent directory name
    return int(k)

def main():
    results = []
    
    # Find all .sam files
    for sam_path in Path('alignment_results/reference').rglob('*.sam'):
        k = extract_k(sam_path)
        score = extract_as_score(sam_path)
        if score is not None:
            results.append({'k': k, 'alignment_score': score})
        else:
            print(f"Failed to extract alignment score from {sam_path}")
    
    # Create and sort DataFrame
    df = pd.DataFrame(results)
    df = df.sort_values('k')
    
    print(df)
    # Optionally save to CSV
    df.to_csv('alignment_scores.csv', index=False)

if __name__ == '__main__':
    main()