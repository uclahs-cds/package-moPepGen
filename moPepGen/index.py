""" moPepGen index """
from __future__ import annotations
from typing import TYPE_CHECKING, Set, List
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

class CanonicalPoolMetadata:
    """ Canonical peptide pool metadata """
    def __init__(self, filename:str, index:int, cleavage_params:CleavageParams):
        """ constructor """
        self.filename = filename
        self.index = index
        self.cleavage_params = cleavage_params

    def jsonfy(self):
        """ jsonfy """
        return {
            'filename': self.filename,
            'index': self.index,
            'cleavage_params': self.cleavage_params.jsonfy()
        }

class IndexMetadata:
    """ Index metadata """
    def __init__(self, version:MetaVersion, canonical_pools:List[CanonicalPoolMetadata],
            source:str):
        """ constructor """
        self.version = version
        self.canonical_pools = canonical_pools
        self.source = source

    def jsonfy(self):
        """ jsonfy """
        return {
            'version': self.version.jsonfy(),
            'canonical_pools': [it.jsonfy() for it in self.canonical_pools],
            'source': self.source
        }

    def register_canonical_pool(self, cleavage_params:CleavageParams):
        """ Register for new canonical pool """
        if self.get_canonical_pool(cleavage_params):
            raise ValueError(
                "Canonical peptide pool already exists with the parameters."
            )
        index = max(it.index for it in self.canonical_pools) + 1 if self.canonical_pools else 1
        pool = CanonicalPoolMetadata(
            filename = f"canonical_peptides_{index:03}.pkl",
            index=index,
            cleavage_params=cleavage_params
        )
        self.canonical_pools.append(pool)
        return pool

    def get_canonical_pool(self, cleavage_params:CleavageParams) -> CanonicalPoolMetadata:
        """ Get canonical peptide pool metadata that match with the given
        cleavage parameters. """
        this = cleavage_params.jsonfy(graph_params=False)
        for pool in self.canonical_pools:
            that = pool.cleavage_params.jsonfy(graph_params=False)
            if this == that:
                return pool
        return None

class IndexDir:
    """ moPepGen index directory """
    def __init__(self, path:Path):
        """ constructor """
        self.path = path
        self.genome_file = self.path/'genome.pkl'
        self.annotation_file = self.path/'annotation.gtf'
        self.proteome_file = self.path/'proteome.pkl'
        self.coding_tx_file = self.path/'coding_transcripts.pkl'
        self.metadata_file = self.path/'metadata.json'
        if self.metadata_file.exists():
            self.metadata = self.load_metadata()
        else:
            self.init_metadata()

    def init_metadata(self):
        """ Initialize metadata """
        self.metadata = IndexMetadata(
            version=MetaVersion(),
            canonical_pools=[],
            source=None
        )

    def load_metadata(self) -> IndexMetadata:
        """ Get metadata of the moPepGen index """
        with open(self.metadata_file, 'rt') as handle:
            data = json.load(handle)
            version = MetaVersion(**data['version'])
            canonical_pools = []
            for it in data['canonical_pools']:
                cleavage_params = CleavageParams(**it['cleavage_params'])
                pool = CanonicalPoolMetadata(
                    filename=it['filename'],
                    index=it['index'],
                    cleavage_params=cleavage_params
                )
                canonical_pools.append(pool)
            source = data['source']
        return IndexMetadata(
            version=version,
            canonical_pools=canonical_pools,
            source=source
        )

    def validate_metadata(self) -> bool:
        """ Valid the index metadata """
        cur_version = MetaVersion()
        if not cur_version.is_valid(self.metadata.version):
            raise err.InvalidIndexError(cur_version, self.metadata.version)
        return True

    def save_metadata(self):
        """ Save metadata to index directory """
        data = self.metadata.jsonfy()
        with open(self.metadata_file, 'wt') as handle:
            json.dump(data, handle, indent=2)

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
        anno = gtf.GenomicAnnotationOnDisk()
        anno.init_handle(self.annotation_file)
        anno.load_index(self.annotation_file, source=self.metadata.source)
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

    def load_canonical_peptides(self, cleavage_params:CleavageParams) -> Set[str]:
        """ Load canonical peptides from index directory """
        pool_data = self.metadata.get_canonical_pool(cleavage_params)
        if not pool_data:
            raise ValueError(
                'No canonical peptide pool match with the cleavage parameters'
            )
        with open(self.path/pool_data.filename, 'rb') as handle:
            return pickle.load(handle)

    def save_canonical_peptides(self, seqs:Set[str], cleavage_params:CleavageParams,
            override:bool=False):
        """ Save canonical peptides to index directory """
        pool_metadata = self.metadata.get_canonical_pool(cleavage_params)
        if not pool_metadata or not override:
            pool_metadata = self.metadata.register_canonical_pool(cleavage_params)
        with open(self.path/pool_metadata.filename, 'wb') as handle:
            pickle.dump(seqs, handle)

    def wipe_canonical_peptides(self):
        """ Remove all canonical peptide pools """
        for pool in self.metadata.canonical_pools:
            os.remove(self.path/pool.filename)
        self.metadata.canonical_pools = []

    def load_coding_tx(self) -> Set[str]:
        """ Load coding transcripts from index directory """
        with open(self.coding_tx_file, 'rb') as handle:
            return pickle.load(handle)

    def save_coding_tx(self, coding_tx:Set[str]):
        """ Save coding transcripts to index directory """
        with open(self.coding_tx_file, 'wb') as handle:
            pickle.dump(coding_tx, handle)
