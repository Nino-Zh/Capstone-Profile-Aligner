import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from tqdm import tqdm 

from align_and_score import (
    getSequences_synthetic,
    gotoh,
    vcf_to_alignment,
    get_mutation_data,
    calculate_total_score,
    seqProfile_modified
)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run alignment analysis with customizable parameters')
    parser.add_argument('--test_cases', nargs='+', type=int, default=[1],
                        help='List of test case numbers to process')
    parser.add_argument('--comparison_groups', nargs='+', type=str,
                        default=["subs_0.01_indel_0.002", "subs_0.05_indel_0.01", "subs_0.09_indel_0.018"],
                        help='Comparison groups to analyze')
    parser.add_argument('--alpha_values', nargs='+', type=float, default=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
                        help='Alpha values to test')
    parser.add_argument('--num_runs', type=int, default=3,
                        help='Number of runs to process per group')
    parser.add_argument('--base_dir', type=str,
                        default="/mnt/c/Users/ninoz_5pamwj8/OneDrive/Desktop/capstone_all/synthetic_data",
                        help='Base directory for data files')
    return parser.parse_args()

def process_group(test_case, group, num_runs, reference_files, base_dir, alpha_values):
    results = []
    parts = group.split('_')
    subs_rate = float(parts[1])
    indel_rate = float(parts[3])
    total_rate = subs_rate + (2 * indel_rate)
    
    # Define genome identifiers consistently
    genomes = ["Drosophila_reference", "Drosophila_synthetic"]
    
    for run_num in range(1, num_runs+1):
        run_path = os.path.join(
            base_dir,
            f"test_{test_case}",
            "input_data",
            group,
            f"run_{run_num}"
        )
        
        # File paths with consistent naming
        vcf_path = os.path.join(run_path, f"variants_{run_num}.vcf")
        alt_fasta_path = os.path.join(run_path, f"mutated_{run_num}.fasta")
        alt_bw_path = os.path.join(run_path, f"tf_profile_{run_num}.bw")
        
        if not all(os.path.exists(f) for f in [vcf_path, alt_fasta_path, alt_bw_path]):
            print(f"Missing files in {run_path}, skipping...")
            continue
            
        try:
            # Generate path mappings for this run
            genomeFnames = {
                "Drosophila_reference": reference_files["fasta"],
                "Drosophila_synthetic": alt_fasta_path
            }
            
            contributionsFnames = {
                "Drosophila_reference": reference_files["bw"],
                "Drosophila_synthetic": alt_bw_path
            }

            # Get true alignment
            alnA_true, alnB_true = vcf_to_alignment(vcf_path, reference_files["fasta"])
            table1_true, table2_true = get_mutation_data(alnA_true, alnB_true)
            
            # Get sequences with explicit parameters
            seqprof, _ = getSequences_synthetic(
                normalize=lambda x: (x/np.sum(x) * x.shape[0]),
                genomes=genomes,
                genomeFnames=genomeFnames,
                contributionsFnames=contributionsFnames
            )
            
            # Process all alpha values
            for alpha in alpha_values:
                alnA_algo, alnB_algo, dist = gotoh(
                    seqprof["Drosophila_reference"],
                    seqprof["Drosophila_synthetic"],
                    u=2,  # gap extension penalty
                    v=5,  # gap open penalty
                    delta=seqProfile_modified(alpha, relative=False, threshold=0)
                )
                
                table1_algo, table2_algo = get_mutation_data(alnA_algo, alnB_algo)
                score = calculate_total_score(table1_algo, table2_algo, table1_true, table2_true)
                
                results.append({
                    "test_case": test_case,
                    "group": group,
                    "total_rate": total_rate,
                    "alpha": alpha,
                    "score": score,
                    "run": run_num
                })
                
        except Exception as e:
            print(f"Error processing {run_path}: {str(e)}")
            
    return results

def generate_plot(df, base_dir, alpha_values, test_cases):
    plt.figure(figsize=(14, 8))
    
    # Calculate statistics using Standard Error
    plot_data = df.groupby(['group', 'total_rate', 'alpha'])['score'].agg(
        mean_score='mean',
        std_err=lambda x: x.std()/np.sqrt(x.count()),
        n='count'
    ).reset_index()

    # Create group labels and positions
    groups = df['group'].unique()
    group_info = {}
    for group in groups:
        parts = group.split('_')
        subs = parts[1]
        indel = parts[3]
        total = df[df['group'] == group]['total_rate'].iloc[0]
        group_info[group] = {
            'label': f"SN: {subs}\nIN: {indel}",
            'total_rate': total
        }

    # Sort groups by total mutation rate
    sorted_groups = sorted(groups, key=lambda x: group_info[x]['total_rate'])
    x_positions = np.arange(len(sorted_groups))
    
    # Create distinct color palette
    colors = plt.cm.gist_ncar(np.linspace(0.1, 0.9, len(alpha_values)))
    
    # Position offsets for alpha values
    n_alpha = len(alpha_values)
    offset = np.linspace(-0.25, 0.25, n_alpha)  # Tight spacing

    # Plot each alpha value
    for idx, alpha in enumerate(sorted(alpha_values)):
        alpha_data = plot_data[plot_data['alpha'] == alpha]
        
        # Get values in group order
        means = []
        errs = []
        for group in sorted_groups:
            group_data = alpha_data[alpha_data['group'] == group]
            means.append(group_data['mean_score'].iloc[0])
            errs.append(group_data['std_err'].iloc[0])
        
        plt.errorbar(
            x=x_positions + offset[idx],
            y=means,
            yerr=errs,
            fmt='o',
            color=colors[idx],
            capsize=4,
            markersize=8,
            label=f'α={alpha}',
            elinewidth=1.5
        )

    # Add group separation lines
    for pos in x_positions[:-1]:
        plt.axvline(pos + 0.5, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    # Configure axes
    plt.xticks(x_positions, [group_info[g]['label'] for g in sorted_groups])
    plt.xlabel("Comparison Groups (Substitution and Indel Rates)")
    plt.ylabel("Average Alignment Score")
    test_case_str = f"Test Case {test_cases[0]}" if len(test_cases) == 1 else f"Test Cases {', '.join(map(str, test_cases))}"
    plt.title(f"{test_case_str}", pad=20, fontsize=14, fontweight='semibold')
    
    # Style adjustments
    plt.grid(True, axis='y', linestyle=':', alpha=0.3)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    
    plt.legend(title="Alpha Values", bbox_to_anchor=(1.1, 1), loc='upper left')
    plt.tight_layout()

    plot_path = os.path.join(base_dir, "mutation_rate_vs_alignment_score.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved to {plot_path}")


def main():
    args = parse_arguments()
    
    reference_files = {
        "fasta": os.path.join(args.base_dir, "FASTA/dm6_reference.fasta"),
        "bw": os.path.join(args.base_dir, "BW/reference_all_neutral.bw")
    }

    all_results = []
    
    for test_case in args.test_cases:
        for group in args.comparison_groups:
            print(f"Processing Test {test_case}, Group {group}...")
            group_results = process_group(test_case=test_case, group=group, num_runs=args.num_runs, reference_files=reference_files,
                                          base_dir=args.base_dir, alpha_values=args.alpha_values)
            all_results.extend(group_results)
    
    df = pd.DataFrame(all_results)
    
    if df.empty:
        raise ValueError("No data collected - check paths and input files!")
    
    generate_plot(df, args.base_dir, args.alpha_values, args.test_cases)

if __name__ == "__main__":
    main()
