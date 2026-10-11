// Multiple scaffold outputs can belong to the same sample or haplotype.
workflow REPORT_SAMPLES {
    take:
    assemblies

    main:
    groups = assemblies
        .map { row -> [sample: [id: row.meta.id, group: row.meta.group]] }
        .unique { row -> row.sample.id }
        .collect()

    emit:
    groups
}
