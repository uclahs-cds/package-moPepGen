""" moPepGen index """
from __future__ import annotations
from typing import TYPE_CHECKING, Set
import os
import shutil
from pathlib import Path
import gzip
import json
import pickle
from moPepGen import gtf, logger, err
from moPepGen.version import MetaVersion
from moPepGen.params import CleavageParams


if TYPE_CHECKING:
    from moPepGen.dna import DNASeqDict
    from moPepGen.gtf import GenomicAnnotationOnDisk
    from moPepGen.aa import AminoAcidSeqDict
    from moPepGen.gtf.GTFPointer import TranscriptPointer, GenePointer

class IndexMetadata:
    """ Index metadata """
    def __init__(self, version:MetaVersion, cleavage_params:CleavageParams,
            source:str):
        """ constructor """
        self.version = version
        self.cleavage_params = cleavage_params
        self.source = source

    def jsonfy(self):
        """ jsonfy """
        return {
            'version': self.version.jsonfy(),
            'cleavage_params': self.cleavage_params.jsonfy(graph_params=False),
            'source': self.source
        }

class IndexDir:
    """ moPepGen index directory """
    def __init__(self, path:Path):
        """ constructor """
        self.path = path
        self.genome_file = self.path/'genome.pkl'
        self.annotation_file = self.path/'annotation.gtf'
        self.proteome_file = self.path/'proteome.pkl'
        self.coding_tx_file = self.path/'coding_transcripts.pkl'
        self.canonical_peptides_file = self.path/'canonical_peptides.pkl'
        self.metadata_file = self.path/'metadata.json'

    def load_genome(self) -> DNASeqDict:
        """ Load genome from index directory """
        with open(self.genome_file, 'rb') as handle:
            return pickle.load(handle)

    def save_genome(self, genome:DNASeqDict):
        """ Save genome to index directory """
        with open(self.genome_file, 'wb') as handle:
            pickle.dump(genome, handle)

    def load_annotation(self) -> GenomicAnnotationOnDisk:
        """ Load genomic annotation from index directory """
        metadata = self.load_metadata()
        anno = gtf.GenomicAnnotationOnDisk()
        anno.init_handle(self.annotation_file)
        anno.load_index(self.annotation_file, source=metadata.source)
        return anno

    def create_gtf_copy(self, file:Path, symlink:bool=True):
        """ Create copy of GTF """
        if file.suffix.lower() == '.gz':
            if symlink:
                symlink = False
                logger(
                    "--gtf-symlink was suppressed because compressed GTF file was received. "
                )
        elif file.suffix.lower() != '.gtf':
            raise ValueError(f"Cannot handle gtf file {file}")

        if symlink:
            os.symlink(file.absolute(), self.annotation_file)
        elif file.suffix.lower() == '.gtf':
            shutil.copy2(file, self.annotation_file)
        else:
            with gzip.open(file, 'rt') as ihandle, open(self.annotation_file, 'wt') as ohandle:
                for line in ihandle:
                    ohandle.write(line)

    def save_annotation(self, file:Path, source:str=None, proteome:AminoAcidSeqDict=None,
            invalid_protein_as_noncoding:bool=True, symlink:bool=True) -> GenomicAnnotationOnDisk:
        """ Save genomic annotation to index directory """
        self.create_gtf_copy(file, symlink=symlink)
        anno = gtf.GenomicAnnotationOnDisk()
        anno.generate_index(self.annotation_file, source)

        if proteome:
            anno.check_protein_coding(proteome, invalid_protein_as_noncoding)

        gene_idx_file, tx_idx_file = anno.get_index_files(self.annotation_file)

        with open(gene_idx_file, 'wt') as handle:
            for gene in anno.genes.keys():
                gene_pointer:GenePointer = anno.genes.get_pointer(gene)
                handle.write(gene_pointer.to_line() + '\n')

        with open(tx_idx_file, 'wt') as handle:
            for tx in anno.transcripts.keys():
                tx_pointer:TranscriptPointer = anno.transcripts.get_pointer(tx)
                handle.write(tx_pointer.to_line() + '\n')
        return anno

    def load_proteome(self) -> AminoAcidSeqDict:
        """ Load proteome from index directory """
        with open(self.proteome_file, 'rb') as handle:
            return pickle.load(handle)

    def save_proteome(self, proteome:AminoAcidSeqDict):
        """ Save proteome to index directory """
        with open(self.proteome_file, 'wb') as handle:
            pickle.dump(proteome, handle)

    def load_canonical_peptides(self) -> Set[str]:
        """ Load canonical peptides from index directory """
        with open(self.canonical_peptides_file, 'rb') as handle:
            return pickle.load(handle)

    def save_canonical_peptides(self, seqs:Set[str]):
        """ Save canonical peptides to index directory """
        with open(self.canonical_peptides_file, 'wb') as handle:
            pickle.dump(seqs, handle)

    def load_coding_tx(self) -> Set[str]:
        """ Load coding transcripts from index directory """
        with open(self.coding_tx_file, 'rb') as handle:
            return pickle.load(handle)

    def save_coding_tx(self, coding_tx:Set[str]):
        """ Save coding transcripts to index directory """
        with open(self.coding_tx_file, 'wb') as handle:
            pickle.dump(coding_tx, handle)

    def load_metadata(self) -> IndexMetadata:
        """ Get metadata of the moPepGen index """
        with open(self.metadata_file, 'rt') as handle:
            data = json.load(handle)
            version = MetaVersion(**data['version'])
            cleavage_params = CleavageParams(**data['cleavage_params'])
            source = data['source']
        return IndexMetadata(
            version=version,
            cleavage_params=cleavage_params,
            source=source
        )

    def validate_metadata(self, cleavage_params:CleavageParams=None) -> bool:
        """ Valid the index metadata """
        metadata = self.load_metadata()
        cur_version = MetaVersion()
        if not cur_version.is_valid(metadata.version):
            raise err.InvalidIndexError(
                cur_version, metadata.version, cleavage_params, metadata.cleavage_params
            )

        if cleavage_params:
            cur_cleavage_json = cleavage_params.jsonfy(graph_params=False)
            index_cleavage_json = metadata.cleavage_params.jsonfy(graph_params=False)
            for k,v in cur_cleavage_json.items():
                if v != index_cleavage_json[k]:
                    raise err.InvalidIndexError(
                        cur_version, metadata.version, cleavage_params, metadata.cleavage_params
                    )
        return True

    def save_metadata(self, metadata:IndexMetadata):
        """ Save metadata to index directory """
        data = metadata.jsonfy()
        with open(self.metadata_file, 'wt') as handle:
            json.dump(data, handle, indent=2)
