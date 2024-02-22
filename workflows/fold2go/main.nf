include { ALPHAFOLD; MSA } from '../../modules/alphafold'
include { PYMOL          } from '../../modules/pymol'
include { COMMS          } from '../../modules/comms'

workflow FOLD2GO {

    Channel
        .fromPath(params.IN)
        .map { fasta -> def (chainA, chainB) = fasta.baseName.minus("${params.ALPHAFOLD.run}_").tokenize("."); [ [ 'A': chainA, 'B': chainB ], fasta ] }
        .set { fasta }

        COMMS(
            fasta.count(),
            file("${params.OUT}/${workflow.runName}", type: 'dir')
        )

        MSA(
            fasta.splitFasta(record: [id: true, sequence: true]).unique{ fasta, record -> record },
            ['uniref90', 'mgnify', 'uniprot', 'bfd']
        )
        
        MSA.out.msa
            .groupTuple(by:0, size:4)
            .combine(fasta)
            .branch { chain, msa, meta, fasta ->
                A: chain == meta.A
                    return [meta, msa]
                B: chain == meta.B
                    return [meta, msa]
            }
            .set { msa }

        ALPHAFOLD(
            msa.A.join(msa.B).join(fasta)
        )
        
        PYMOL(
            ALPHAFOLD.out.prediction.combine(fasta, by:0)
        )

        PYMOL.out.metrics.collectFile(name: "template_indep_info.tsv", storeDir: "${params.OUT}/${workflow.runName}", keepHeader: true)
}
