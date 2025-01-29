import pandas as pd
from pathlib import Path
from subprocess import run, PIPE
import re

Path('alignment_results').mkdir(exist_ok=True)

for tech in ['ONT', 'PB']:
    df = []
    for concat in [True]:
        ref = f"ecoli/ecoli_{tech}{'_concat' if concat else ''}.fasta"
        samples_dir = Path(f"samples/{tech}/{'concat' if concat else 'non_concat'}")
        
        if not Path(ref).exists() or not samples_dir.exists():
            continue  # Skip missing paths

        # Set presets
        if tech == 'ONT':
            preset = '-ax map-ont'
        else:  # ONT or PB
            preset = f'-ax map-pb'

        # Process samples with silent output
        for fasta in samples_dir.glob('*/*.fasta'):
            k = fasta.parent.name
            output_dir = Path(f"sam_files/{tech}/{'concat' if concat else 'non_concat'}/{k}")
            output_dir.mkdir(parents=True, exist_ok=True)
            output_file = output_dir / f"{fasta.stem}.sam"
                            
            command = f'./minimap2/minimap2 -t 128 {preset} {ref} {fasta} > {output_file} 2>/dev/null'
            print(command)
            if not output_file.exists():
                result = run(
                    command,
                    shell=True, stdout=PIPE, text=True
                )
            else:
                print(f"Skipping existing file {output_file}")
                
            grep_cmd = f"grep 'AS:i:' {output_file} | awk -F'AS:i:' '{{print $2}}' | awk '{{print $1}}'"
            grep_result = run(grep_cmd, shell=True, stdout=PIPE, text=True)
            scores = [int(x) for x in grep_result.stdout.splitlines() if x]
            max_score = max(scores) if scores else None
            print(f"Max score is: {max_score}")

            df.append({
                'k': int(k),
                'score': max_score,
                'concatenated': concat
            })

    
    # Save results
    if df:
        pd.DataFrame(df).to_csv(f'alignment_results/{tech.lower()}.csv', index=False)

ref = "ecoli/ref.fasta"
samples_dir = Path("samples/reference/")
df_ref = []

if Path(ref).exists() and samples_dir.exists():
    # Use assembly-to-assembly preset for reference comparison
    preset = '-ax asm5'  
    
    for fasta in samples_dir.glob("*/*.fasta"):
        k = fasta.parent.name
        output_dir = Path(f"sam_files/reference/{k}")
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{fasta.stem}.sam"
        
        # Run minimap2 command
        command = f'./minimap2/minimap2 -t 128 {preset} {ref} {fasta} > {output_file} 2>/dev/null'
        print(command)
        
        if not output_file.exists():
            result = run(command, shell=True, stdout=PIPE, text=True)
        
        # Extract alignment scores
        grep_cmd = f"grep 'AS:i:' {output_file} | awk -F'AS:i:' '{{print $2}}' | awk '{{print $1}}'"
        grep_result = run(grep_cmd, shell=True, stdout=PIPE, text=True)
        scores = [int(x) for x in grep_result.stdout.splitlines() if x]
        max_score = max(scores) if scores else None
        
        df_ref.append({
            'k': int(k),
            'score': max_score,
            'sample': fasta.stem
        })

# Save reference results
if df_ref:
    pd.DataFrame(df_ref).to_csv('alignment_results/reference.csv', index=False)
    print("Saved reference alignment results")
else:
    print("No reference alignment results to save")
