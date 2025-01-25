import os
from pathlib import Path

# Reference genome path
REF_GENOME = "ecoli/ref.fasta"
# Input and output root directories
INPUT_DIR = "samples/reference"
OUTPUT_DIR = "alignment_results/reference"

def main():
    # Create output directory if it doesn't exist
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    
    # Find all .fasta files recursively
    for fasta_path in Path(INPUT_DIR).rglob("*.fasta"):
        # Get relative path from INPUT_DIR
        rel_path = fasta_path.relative_to(INPUT_DIR)
        # Construct output path with .sam extension
        output_path = Path(OUTPUT_DIR) / rel_path.with_suffix('.sam')
        
        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if output_path.exists():
            print(f"Skipping: {fasta_path} -> {output_path}")
            continue
        
        # Run minimap2
        cmd = [
            "./minimap2/minimap2",
            "-ax asm5",
            REF_GENOME,
            str(fasta_path),
            ">",
            str(output_path)
        ]
        
        print(f"Processing: {fasta_path} -> {output_path}")
        os.system(" ".join(cmd) + " 2>/dev/null")

if __name__ == "__main__":
    main()