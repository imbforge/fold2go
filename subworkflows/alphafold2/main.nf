nextflow.enable.types = true

include { MSA; INFERENCE } from '../../modules/alphafold2'

workflow ALPHAFOLD2 {

    take:
        input: Channel<Record>

    main:
        List databases = ['uniref90', 'mgnify', 'bfd', 'uniprot']        

        sequences =
            input
            .map { job ->
                def letters = ('A'..'H')
                record(
                    chains: job.input.splitFasta ( record: [ id: true ] ).withIndex().collectEntries { it, idx -> [(letters[idx]): (it.id)] },
                    fasta:  job.input,
                    model:  job.model
                )
            }

        jobdef =
            sequences
            .combine (
                channel.fromList(databases)
            )
            .flatMap { complex, db ->
                complex.fasta
                .splitFasta ( file: true )
                .collect { chunk ->
                    tuple( complex.chains, chunk.splitFasta ( record: [ id: true] ).id.pop(), chunk, complex.model, db )
                }
            }
            // FIXME
            // this whole tuple detour seems wrong but the typechecker complains if the record is constructed in that .flatMap() above :(
            // same goes for the List unpacking below, not sure why 'map { chains, id, chunk, model, db -> ... }' isn't allowed
            .map { R ->
                def (chains, id, chunk, model, db) = R
                record(
                    chains: chains,
                    id:     id,
                    fasta:  chunk,
                    db:     db,
                    model:  model
                )
            }
            .unique { it -> it.id + it.db }
        
        msa = MSA(jobdef)
        
        prediction = INFERENCE(
            msa
            .flatMap { record ->
                [ tuple(record.chains, databases.size(), record.msa) ]
            }
            .groupBy()
            .map { chains, msas -> 
                record(
                    chains: chains,
                    A: msas.findAll { it -> it.parent.name == (chains['A'] as String) }.toSet(),
                    B: msas.findAll { it -> it.parent.name == (chains['B'] as String) }.toSet(),
                    C: msas.findAll { it -> it.parent.name == (chains['C'] as String) }.toSet(),
                    D: msas.findAll { it -> it.parent.name == (chains['D'] as String) }.toSet(),
                    E: msas.findAll { it -> it.parent.name == (chains['E'] as String) }.toSet(),
                    F: msas.findAll { it -> it.parent.name == (chains['F'] as String) }.toSet(),
                    G: msas.findAll { it -> it.parent.name == (chains['G'] as String) }.toSet(),
                    H: msas.findAll { it -> it.parent.name == (chains['H'] as String) }.toSet()
                )
            }
            .join(sequences, by: 'chains')
        )

    emit:
        msa: Channel<Record> = msa
        prediction: Channel<Record> = prediction
}