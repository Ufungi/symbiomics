//
// Sequence-feature predictors for `-entry functional`: secretion signals,
// subcellular localisation, transmembrane topology, effector-likeness,
// peptidases. Each process wraps ONE command and always emits its output file;
// the file is EMPTY when the tool is not installed (on_missing_tool=skip) or
// the module was not requested, so the parser treats "not run" and "ran, found
// nothing" differently and the pipeline never dies on a licence-gated tool.
//
// Licensing: SignalP, TargetP, DeepLoc, DeepTMHMM and Phobius are DTU/SBC
// academic-licence software and cannot be bundled or auto-downloaded. Install
// them yourself (see docs/functional_annotation.md) and either put the command
// on PATH or set --<tool>_cmd to an absolute path. Command lines below follow
// each tool's own documented CLI.
//

def skip_guard(cmd, out, pname, fail) {
    return """
    if ! command -v ${cmd.tokenize()[0]} >/dev/null 2>&1; then
        echo "[symbiomics] '${cmd.tokenize()[0]}' not found -- ${pname} SKIPPED (empty output)" >&2
        ${fail ? 'exit 127' : ": > ${out}; echo '\"${pname}\": {tool: skipped}' > versions.yml; exit 0"}
    fi
    """
}

def opt_flag(flag, f) {
    return f.size() > 0 ? "${flag} '${f}'" : ''
}

process SIGNALP6 {
    label 'process_medium'
    tag "${proteome.name}"
    input:  path proteome
    output: path 'signalp6.tsv', emit: tsv
            path 'versions.yml', emit: versions
    script:
    """
    ${skip_guard(params.signalp6_cmd, 'signalp6.tsv', task.process, params.on_missing_tool == 'fail')}
    ${params.signalp6_cmd} --fastafile '${proteome}' --organism ${params.signalp_organism} \\
        --output_dir sp6 --format none --mode fast --write_procs ${task.cpus}
    cp sp6/prediction_results.txt signalp6.tsv
    echo '"${task.process}": {signalp: 6}' > versions.yml
    """
    stub:
    """
    printf '# ID\\tPrediction\\tOTHER\\tSP(Sec/SPI)\\tCS Position\\n' > signalp6.tsv
    echo '"${task.process}": {signalp: stub}' > versions.yml
    """
}

process TARGETP2 {
    label 'process_medium'
    tag "${proteome.name}"
    input:  path proteome
    output: path 'targetp2.tsv', emit: tsv
            path 'versions.yml', emit: versions
    script:
    // -org pl (plant: also cTP) or non-pl (fungi/animal: SP, mTP only)
    def org = params.taxon == 'plant' ? 'pl' : 'non-pl'
    """
    ${skip_guard(params.targetp_cmd, 'targetp2.tsv', task.process, params.on_missing_tool == 'fail')}
    ${params.targetp_cmd} -fasta '${proteome}' -org ${org} -format short -prefix tp2
    cp tp2_summary.targetp2 targetp2.tsv
    echo '"${task.process}": {targetp: 2.0}' > versions.yml
    """
    stub:
    """
    printf '# ID\\tPrediction\\tnoTP\\tSP\\tmTP\\tCS Position\\n' > targetp2.tsv
    echo '"${task.process}": {targetp: stub}' > versions.yml
    """
}

process DEEPLOC2 {
    label 'process_medium'
    tag "${proteome.name}"
    input:  path proteome
    output: path 'deeploc2.csv', emit: csv
            path 'versions.yml', emit: versions
    script:
    """
    ${skip_guard(params.deeploc2_cmd, 'deeploc2.csv', task.process, params.on_missing_tool == 'fail')}
    ${params.deeploc2_cmd} -f '${proteome}' -o dl2 -m ${params.deeploc2_model}
    cp dl2/results_*.csv deeploc2.csv
    echo '"${task.process}": {deeploc: 2}' > versions.yml
    """
    stub:
    """
    printf 'Protein_ID,Localizations,Signals\\n' > deeploc2.csv
    echo '"${task.process}": {deeploc: stub}' > versions.yml
    """
}

process DEEPTMHMM {
    label 'process_medium'
    tag "${proteome.name}"
    input:  path proteome
    output: path 'deeptmhmm.3line', emit: topo
            path 'versions.yml'    , emit: versions
    script:
    // DeepTMHMM is distributed through BioLib (pybiolib); a local install that
    // exposes `deeptmhmm`-compatible CLI can be substituted via --deeptmhmm_cmd.
    """
    ${skip_guard(params.deeptmhmm_cmd, 'deeptmhmm.3line', task.process, params.on_missing_tool == 'fail')}
    ${params.deeptmhmm_cmd} --fasta '${proteome}'
    cp biolib_results/predicted_topologies.3line deeptmhmm.3line
    echo '"${task.process}": {deeptmhmm: 1.0}' > versions.yml
    """
    stub:
    """
    : > deeptmhmm.3line
    echo '"${task.process}": {deeptmhmm: stub}' > versions.yml
    """
}

process PHOBIUS {
    label 'process_single'
    tag "${proteome.name}"
    input:  path proteome
    output: path 'phobius.tsv'  , emit: tsv
            path 'versions.yml' , emit: versions
    script:
    """
    ${skip_guard(params.phobius_cmd, 'phobius.tsv', task.process, params.on_missing_tool == 'fail')}
    ${params.phobius_cmd} -short '${proteome}' > phobius.tsv
    echo '"${task.process}": {phobius: 1.01}' > versions.yml
    """
    stub:
    """
    printf 'SEQENCE ID\\tTM\\tSP\\tPREDICTION\\n' > phobius.tsv
    echo '"${task.process}": {phobius: stub}' > versions.yml
    """
}

process EFFECTORP {
    label 'process_medium'
    tag "${proteome.name}"
    input:  path proteome
    output: path 'effectorp.tsv', emit: tsv
            path 'versions.yml' , emit: versions
    script:
    // EffectorP is meant for the secretome, not the whole proteome; the
    // parser only reports it for proteins that are also called secreted.
    """
    ${skip_guard(params.effectorp_cmd, 'effectorp.tsv', task.process, params.on_missing_tool == 'fail')}
    ${params.effectorp_cmd} -i '${proteome}' -o effectorp.tsv
    echo '"${task.process}": {effectorp: 3}' > versions.yml
    """
    stub:
    """
    printf '# Identifier\\tPrediction\\n' > effectorp.tsv
    echo '"${task.process}": {effectorp: stub}' > versions.yml
    """
}

process MEROPS_SCAN {
    label 'process_medium'
    tag "${proteome.name}"
    input:  path proteome
            path merops_dmnd
    output: path 'merops.tsv'   , emit: tsv
            path 'versions.yml' , emit: versions
    script:
    """
    diamond blastp -q '${proteome}' -d '${merops_dmnd}' -o merops.tsv \\
        --outfmt 6 qseqid sseqid pident evalue bitscore stitle \\
        --evalue 1e-5 --max-target-seqs 1 --threads ${task.cpus}
    echo '"${task.process}": {diamond: merops_scan}' > versions.yml
    """
    stub:
    """
    : > merops.tsv
    echo '"${task.process}": {diamond: stub}' > versions.yml
    """
}

process CONSOLIDATE_PREDICTORS {
    label 'process_single'
    tag "${proteome.name}"
    input:  path proteome
            path signalp6
            path targetp2
            path deeploc2
            path deeptmhmm
            path phobius
            path effectorp
            path merops
    output: path 'localization.tsv', emit: tsv
            path 'versions.yml'    , emit: versions
    script:
    // Same sentinel convention as ASSIGN_PRODUCTS: presence is tested by size.
    """
    parse_predictors.py --proteome '${proteome}' \\
        ${opt_flag('--signalp6', signalp6)} ${opt_flag('--targetp2', targetp2)} ${opt_flag('--deeploc2', deeploc2)} \\
        ${opt_flag('--deeptmhmm', deeptmhmm)} ${opt_flag('--phobius', phobius)} ${opt_flag('--effectorp', effectorp)} \\
        ${opt_flag('--merops', merops)} --secretome-min-tools ${params.secretome_min_tools} --out localization.tsv
    echo '"${task.process}": {python: '\$(python3 --version | sed 's/Python //')'}' > versions.yml
    """
    stub:
    """
    printf 'protein_id\\tsecreted\\tn_sp_tools\\tn_tm\\tlocalization\\teffector\\tmerops\\n' > localization.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
