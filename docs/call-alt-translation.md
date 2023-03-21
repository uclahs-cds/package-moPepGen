{% set command = 'callAltTranslation' %}
# {{ command }}

::: moPepGen.cli.call_alt_translation
	handler: python
    selection:
      members: false
    rendering:
      show_root_heading: false
      show_source: false

{% include 'partials/_caution_on_reference_version.md' %}

## Usage

```
{{ get_arg_usage(command) }}
```

## Arguments

{% with actions=get_arg_data(command) %}
{% include 'partials/_command_usage.md' %}
{% endwith %}

## Alternative translation

Alternative translation is when a different peptide is generated from the same transcript without any change of the genetic code from the genome.

### Selenocysteine Termination

In eukaryotes, the UGA on some mRNAs can be decoded into selenocysteine instead of being recognized as a stop codon, and those proteins are called selenoproteins. However the decoding of UGA is regulated by complex signals including mRNA and sec-tRNA abundance, which could result two isoforms: one with UGA read through and one being truncated. Selenocysteine termination is used to represent the later situation. Selenocysteine terminations are not written into any GVF file but they are represented in the format of `SECT-<pos>` where `pos` is the position of the selenocysteine UGA being recognized as a stop codon in the **gene**.

### Tryptophan > Phenylalanine Codon Reassignment

Tryptophan > Phenylalanine substitutants, described in [Patasker, et al.](https://pubmed.ncbi.nlm.nih.gov/35264796/), happens when cellular tryptophan is depleted and phenylalanine is reassigned to tryptophan codons to have protein synthesis continue. The process largely exists in tumor cells. Similar to selenocysteine termination, W > F substitutants are also not written into GVFs, but is represented in the format of `W2F-<pos>`. Noted that the `pos` is at peptide coordinate (*i.e.,* zeroed at the beginning of the peptide).
