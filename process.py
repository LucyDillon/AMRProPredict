#!/usr/bin/env python3
import sys
import argparse
import os
import shutil
import subprocess
import yaml
from pathlib import Path
import tempfile
import tarfile
import gc

def check_singularity_available():
    """Check if Singularity is available"""
    try:
        subprocess.run(['singularity', '--version'], 
                       capture_output=True, text=True, check=True)
        return True
    except Exception:
        return False

def safe_copy(src, dst):
    """Copy file safely"""
    try:
        Path(dst).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        return Path(dst).exists()
    except Exception as e:
        print(f"Copy failed: {e}", file=sys.stderr)
        return False

def find_file(tool_dir, filename):
    """Find a file - returns path or None"""
    tool_path = Path(tool_dir)
    if (tool_path / filename).exists():
        return str(tool_path / filename)
    if (tool_path / 'data' / filename).exists():
        return str(tool_path / 'data' / filename)
    return None

def create_minimal_config(work_dir, input_basename, tool_dir, model_type):
    """Create minimal config and copy required files to working directory"""
    tool_path = Path(tool_dir)

    # Find required files
    diamond_db = find_file(tool_dir, 'eggnog_proteins.dmnd') or str(tool_path / 'eggnog_proteins.dmnd')
    eggnog_db = find_file(tool_dir, 'eggnog.db') or str(tool_path / 'eggnog.db')

    # Copy essential scripts to working directory
    scripts_to_copy = {
        'eggnog_output_to_arff.py': None,
        'eggnog2cnnPredictions.py': None,
        'genefamilies.txt': None
    }
    for script_name in scripts_to_copy:
        src = find_file(tool_dir, script_name)
        if src:
            dst = work_dir / script_name
            if safe_copy(src, dst):
                scripts_to_copy[script_name] = str(dst)
                print(f"Copied {script_name} to working directory", file=sys.stderr)
            else:
                scripts_to_copy[script_name] = src
        else:
            scripts_to_copy[script_name] = str(tool_path / script_name)

    # Copy model files to working directory
    model_files_dir = tool_path / 'data'
    model_files_copied = {}
    
    # Define which files to copy based on model type
    if model_files_dir.exists():
        # Copy .dot files for decision tree models (small files - avoid segfault with long paths)
        if model_type == 'decision_tree' or model_type == 'both_models':
            for dot_file in model_files_dir.glob('*.dot'):
                dst = work_dir / dot_file.name
                if safe_copy(dot_file, dst):
                    model_files_copied[dot_file.name] = str(dst)
                    print(f"Copied model file: {dot_file.name}", file=sys.stderr)
                else:
                    model_files_copied[dot_file.name] = str(dot_file)
                    print(f"Warning: Failed to copy {dot_file.name}", file=sys.stderr)
        
        # Copy CNN model files if needed (also small files)
        if model_type == 'cnn' or model_type == 'both_models':
            for ext in ['*.h5', '*.pkl', '*.json']:
                for model_file in model_files_dir.glob(ext):
                    dst = work_dir / model_file.name
                    if safe_copy(model_file, dst):
                        model_files_copied[model_file.name] = str(dst)
                        print(f"Copied model file: {model_file.name}", file=sys.stderr)
                    else:
                        model_files_copied[model_file.name] = str(model_file)
                        print(f"Warning: Failed to copy {model_file.name}", file=sys.stderr)
    
    # Copy Singularity containers
    singularity_dir = work_dir / "singularity"
    singularity_dir.mkdir(exist_ok=True)

    container_sources = {
        'container': tool_path / 'amrmlpipeline_DT_arm64.sif',
        'ml_container': tool_path / 'amrmlpipeline_CNN.sif'
    }
    
    # Also check in parent directory
    if not container_sources['container'].exists():
        container_sources['container'] = Path('/__YOUR__PATH___/amrmlpipeline_DT_arm64.sif')
    if not container_sources['ml_container'].exists():
        container_sources['ml_container'] = Path('/__YOUR__PATH___/amrmlpipeline_CNN.sif')

    singularity_copies = {}
    for name, src in container_sources.items():
        if src.exists():
            dest = singularity_dir / src.name
            if safe_copy(src, dest):
                singularity_copies[name] = str(dest)
                print(f"Copied Singularity container: {src.name}", file=sys.stderr)
            else:
                singularity_copies[name] = str(src)
                print(f"Warning: Failed to copy {src.name}, using original path", file=sys.stderr)
        else:
            singularity_copies[name] = str(src)
            print(f"Warning: Missing container: {src}", file=sys.stderr)

    # Build config - use working directory for .dot files, original location for large databases
    config = {
        'genome': input_basename,
        'diamond_db': diamond_db,
        'eggnog_db': eggnog_db,
        'eggnog_data_dir': str(tool_path),
        'eggnog_2_arff': scripts_to_copy['eggnog_output_to_arff.py'],
        'eggnog_2_cnn': scripts_to_copy['eggnog2cnnPredictions.py'],
        'gene_families': scripts_to_copy['genefamilies.txt'],
        'model_files_dir': str(work_dir),  # For .dot files (copied locally)
        'original_model_files_dir': str(model_files_dir),  # For large databases (not copied)
        'tool_dir': str(tool_path),
        'singularity': {
            'enabled': check_singularity_available(),
            'container': singularity_copies['container'],
            'ml_container': singularity_copies['ml_container']
        }
    }

    config_path = work_dir / 'config.yaml'
    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    print("Minimal config created", file=sys.stderr)
    print(f"Model files dir (for .dot): {config['model_files_dir']}", file=sys.stderr)
    print(f"Original model files dir (for .sbt.zip): {config['original_model_files_dir']}", file=sys.stderr)
    return config_path, singularity_dir


def run_snakemake_workflow(input_file, output_file, model_type, cores=1, tool_dir=None):
    """Run Snakemake workflow with Singularity container support"""
    if not tool_dir:
        tool_dir = Path(__file__).parent.absolute()
    else:
        tool_dir = Path(tool_dir)

    snakefiles = {
        'decision_tree': 'decision_tree.smk',
        'cnn': 'cnn.smk',
        'both_models': 'both_models.smk'
    }
    snakefile_path = tool_dir / snakefiles[model_type]
    if not snakefile_path.exists():
        raise FileNotFoundError(f"Snakefile not found: {snakefile_path}")

    temp_dir = tempfile.mkdtemp(prefix='snakemake_')
    try:
        work_dir = Path(temp_dir)
        print(f"Working in: {work_dir}", file=sys.stderr)

        # Copy input
        input_path = Path(input_file)
        input_basename = input_path.stem
        input_copy = work_dir / f"{input_basename}.fna"
        if not safe_copy(input_file, input_copy):
            raise IOError("Failed to copy input file")

        # Create config and copy required files (including .dot files)
        config_path, singularity_dir = create_minimal_config(work_dir, input_basename, tool_dir, model_type)
        gc.collect()

        cmd = [
            'snakemake',
                '--snakefile', str(snakefile_path),
                '--configfile', str(config_path),
                '--cores', str(cores),
                '--directory', str(work_dir),
                '--printshellcmds',
                '--rerun-incomplete'
        ]

        # Add Singularity support
        if check_singularity_available():
            bind_paths = f"{work_dir}:{work_dir},{tool_dir}:{tool_dir},{singularity_dir}:{singularity_dir}"
            cmd.extend([
                '--use-singularity',
                '--singularity-args',
                f'--cleanenv --contain --bind {bind_paths}'
            ])

        # Setup environment
        env = os.environ.copy()
        env['TMPDIR'] = str(work_dir / 'tmp')
        env['SINGULARITY_TMPDIR'] = str(work_dir / 'tmp')
        env['SINGULARITY_CACHEDIR'] = str(work_dir / 'singularity_cache')
        os.makedirs(env['TMPDIR'], exist_ok=True)
        os.makedirs(env['SINGULARITY_CACHEDIR'], exist_ok=True)

        print("Running Snakemake...", file=sys.stderr)
        with open(work_dir / 'snakemake.log', 'w') as log_file:
            subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, check=True, env=env)

        print("Snakemake completed", file=sys.stderr)

        # Package output
        expected_output = work_dir / f"{input_basename}_results.tar.gz"
        if not expected_output.exists():
            result_files = [
                f for f in work_dir.iterdir()
                if f.is_file() and f.suffix in ['.arff', '.txt', '.tsv', '.annotations']
            ]
            if result_files:
                with tarfile.open(expected_output, "w:gz") as tar:
                    for f in result_files:
                        tar.add(f, arcname=f.name)

        if expected_output.exists():
            shutil.copy2(expected_output, output_file)
            print(f"Output created: {output_file}", file=sys.stderr)
            return True
        else:
            print("Snakemake log:", file=sys.stderr)
            print((work_dir / 'snakemake.log').read_text(), file=sys.stderr)
            raise FileNotFoundError("No output generated")

    except subprocess.CalledProcessError:
        print("Snakemake failed", file=sys.stderr)
        log_file = work_dir / 'snakemake.log'
        if log_file.exists():
            print("--- Snakemake log start ---", file=sys.stderr)
            with open(log_file, 'r', errors='ignore') as lf:
                for line in lf:
                    sys.stderr.write(line)
            print("--- Snakemake log end ---", file=sys.stderr)
        raise



def main():
    parser = argparse.ArgumentParser(description='Run ML analysis')
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--model', required=True, 
                        choices=['decision_tree', 'cnn', 'both_models'])
    parser.add_argument('--cores', type=int, default=1)
    parser.add_argument('--tool-dir', required=True)
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print("Error: Input file not found", file=sys.stderr)
        sys.exit(1)

    try:
        run_snakemake_workflow(args.input, args.output, args.model, args.cores, args.tool_dir)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
