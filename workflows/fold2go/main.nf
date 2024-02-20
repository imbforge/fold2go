include { ALPHAFOLD; MSA } from '../../modules/alphafold'
include { PYMOL          } from '../../modules/pymol'
include { COMMS          } from '../../modules/comms'

workflow FOLD2GO {

    Channel
        .fromPath(params.IN)
        .map { fasta -> def (chainA, chainB) = fasta.baseName.minus("${params.ALPHAFOLD.run}_").tokenize("."); [ ['A': chainA, 'B': chainB], fasta ] }
        .tap { fasta }
        .splitFasta(record: [id: true, sequence: true])
        .unique{ chains, fasta -> fasta.sequence }
        .collectFile(cache: true, storeDir: "${params.OUT}/chains"){ chains, fasta -> ["${fasta.id}.fasta", [">chain_${chains.find{ it.value == fasta.id }.key}",fasta.sequence].join("\n")] }
        .set { chains }

        MSA(
            chains,
            ['uniref90', 'mgnify', 'uniprot', 'bfd']
        )

        MSA.out.msa
            .groupTuple(by:0, size: 4)
            .combine(fasta)
            .branch { id, fasta, chains, msa ->
                A: id == chains.A
                    return [chains*.value.join("."), fasta, msa]
                B: id == chains.B
                    return [chains*.value.join("."), fasta, msa]
            }
            .set{ msa }
 
        ALPHAFOLD(
            msa.A.combine(msa.B, by:[0,2])
        )
 
        PYMOL(
            ALPHAFOLD.out.prediction
        )

        COMMS(
            chains.count(),
            file("${params.OUT}/${workflow.runName}", type: 'dir')
        )

}
