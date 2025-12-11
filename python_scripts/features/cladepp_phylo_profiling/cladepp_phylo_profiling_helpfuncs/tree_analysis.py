import os
import subprocess


def analyze_tree_clades(
    tree_path,
    mgc_candidates_dir,
    mapping_file,
    compute_gain_loss_coevolution,
    max_workers=32,
    use_processes=False
):
    if not os.path.exists(mgc_candidates_dir):
        raise ValueError(f"MGC directory does not exist: {mgc_candidates_dir}")

    print(f"Scanning {mgc_candidates_dir}...")

    csv_paths = []
    for sub in sorted(os.listdir(mgc_candidates_dir)):
        subdir = os.path.join(mgc_candidates_dir, sub)
        if not os.path.isdir(subdir):
            continue

        csv = os.path.join(subdir, "comparison_results.csv")
        if os.path.exists(csv):
            csv_paths.append(csv)
            print(f"  ✔ Found: {csv}")
        else:
            print(f"  Skipped: {sub} (no comparison_results.csv)")

    if not csv_paths:
        raise RuntimeError("No comparison_results.csv found!")

    print(f"\nTotal candidates: {len(csv_paths)}")

    # ----------------------------------------------------------------------
    # Write array_input.txt
    # ----------------------------------------------------------------------
    array_file = os.path.join(mgc_candidates_dir, "array_input.txt")
    with open(array_file, "w", newline="\n") as f:
        for p in csv_paths:
            f.write(p + "\n")

    print(f"✔ array_input.txt written: {array_file}")

    # ----------------------------------------------------------------------
    # Create array sbatch files
    # ----------------------------------------------------------------------
    logs_dir = os.path.join(mgc_candidates_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)

    chunk_size = 500            # SAFE: each array < 1000
    num_tasks = len(csv_paths)

    chunks = [(i, min(i+chunk_size-1, num_tasks))
              for i in range(1, num_tasks+1, chunk_size)]
    
    # Skip the first chunk
    if chunks:
        chunks = chunks[1:]

    prev_jobid = None

    for idx, (start, end) in enumerate(chunks, 1):

        array_script = os.path.join(mgc_candidates_dir,
                                    f"compare_array_{idx}.sbatch")

        chunk_len = end - start + 1

        # --------------------------
        # CORRECTED ARRAY TEMPLATE
        # --------------------------
        script_content = f"""#!/bin/bash
#SBATCH -p itaym-pool
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task={max_workers}
#SBATCH --mem=32G
#SBATCH --array=1-{chunk_len}%15
#SBATCH -J treeAnalyzer_{idx}
#SBATCH -o {logs_dir}/%A_%a.out
#SBATCH -e {logs_dir}/%A_%a.err

# Local index inside this chunk (1..chunk_len)
LOCAL_ID=$SLURM_ARRAY_TASK_ID

# Global index in the full list
GLOBAL_ID=$(( LOCAL_ID + {start} - 1 ))

comparison_csv=$(sed -n "${{GLOBAL_ID}}p" {array_file})
mgc_dir=$(dirname "$comparison_csv")

bash /groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_comparison_csv_tree_analyzer.sh \\
    "{tree_path}" \\
    "$comparison_csv" \\
    "$mgc_dir" \\
    "{mapping_file}" \\
    "{str(compute_gain_loss_coevolution).lower()}" \\
    "{max_workers}" \\
    "{str(use_processes).lower()}"
"""

        # Write the sbatch file
        with open(array_script, "w", newline="\n") as f:
            f.write(script_content)

        print(f"✔ Created: {array_script}")

        # --------------------------
        # SUBMIT WITH DEPENDENCY
        # --------------------------
        cmd = ["sbatch"]
        if prev_jobid is not None:
            cmd.append(f"--dependency=afterany:{prev_jobid}")
        cmd.append(array_script)

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            out = result.stdout.strip()
            print(f"  → Submitted: {out}")

            if "Submitted batch job" in out:
                prev_jobid = out.split()[-1]

        except Exception as e:
            print(f"❌ Failed submitting {array_script}: {e}")
